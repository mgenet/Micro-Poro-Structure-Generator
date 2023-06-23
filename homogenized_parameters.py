def homogenization(fname, Em, num, mesh_corners=None, vol=None):

    # from __future__ import print_func/tion
    import dolfin as df
    import math
    import numpy as np
    import meshio
    import matplotlib.pyplot as plt
    ################################################################################
    mesh = df.Mesh()
    df.XDMFFile(fname+"-mesh.xdmf").read(mesh)
    # mesh = df.Mesh(fname + ".xml")
    # mesh = meshio.read(fname+"-mesh.vtk")
    # mesh.points = mesh.points[:, :2]
    # meshio.write(fname+"-mesh.xdmf", mesh)

    # mesh = df.Mesh()
    # df.XDMFFile(fname+"-mesh.xdmf").read(mesh)
    # subdomains = df.MeshFunction("size_t", mesh, fname + "_physical_region.xml")
    # facet = df.MeshFunction("size_t", mesh, fname + "_facet_region.xml")

    # df.XDMFFile("facets.xdmf").write(facet)

    # xdmf_file = df.XDMFFile("subdomains.xdmf")
    # xdmf_file.write(subdomains)
    # xdmf_file.close()

    if mesh_corners is None:
        coord = mesh.coordinates()
        xmax = max(coord[:,0]); xmin = min(coord[:,0])
        ymax = max(coord[:,1]); ymin = min(coord[:,1])
        vol = (ymax - ymin)*(xmax - xmin)
        
        mesh_corners = np.array([[xmin, ymin],
                            [xmax, ymin],
                            [xmax, ymax],
                            [xmin, ymax]])


    # class used to define the periodic boundary map
    class PeriodicBoundary(df.SubDomain):
        def __init__(self, mesh_corners, tolerance=1e-4):# tolerance=df.DOLFIN_EPS):
            """ mesh_corners stores the coordinates of the 4 unit cell corners"""
            df.SubDomain.__init__(self, tolerance)
            self.tol = tolerance
            self.vv = mesh_corners
            self.a1 = self.vv[1,:]-self.vv[0,:] # first vector generating periodicity
            self.a2 = self.vv[3,:]-self.vv[0,:] # second vector generating periodicity
            # check if UC mesh_corners form indeed a parallelogram
            assert np.linalg.norm(self.vv[2, :]-self.vv[3, :] - self.a1) <= self.tol
            assert np.linalg.norm(self.vv[2, :]-self.vv[1, :] - self.a2) <= self.tol

        def inside(self, x, on_boundary):
            # return True if on left or bottom boundary AND NOT on one of the
            # bottom-right or top-left mesh_corners
            return bool((df.near(x[0], self.vv[0,0] + x[1]*self.a2[0]/self.vv[3,1], self.tol) or
                        df.near(x[1], self.vv[0,1] + x[0]*self.a1[1]/self.vv[1,0], self.tol)) and
                        (not ((df.near(x[0], self.vv[1,0], self.tol) and df.near(x[1], self.vv[1,1], self.tol)) or
                        (df.near(x[0], self.vv[3,0], self.tol) and df.near(x[1], self.vv[3,1], self.tol)))) and on_boundary)

        def map(self, x, y):
            if df.near(x[0], self.vv[2,0], self.tol) and df.near(x[1], self.vv[2,1], self.tol): # if on top-right corner
                y[0] = x[0] - (self.a1[0]+self.a2[0])
                y[1] = x[1] - (self.a1[1]+self.a2[1])
            elif df.near(x[0], self.vv[1,0] + x[1]*self.a2[0]/self.vv[2,1], self.tol): # if on right boundary
                y[0] = x[0] - self.a1[0]
                y[1] = x[1] - self.a1[1]
            else:   # should be on top boundary
                y[0] = x[0] - self.a2[0]
                y[1] = x[1] - self.a2[1]

    material_parameters = [(Em, num)]
    nphases = len(material_parameters)

    def eps(v):
        return df.sym(df.grad(v))

    def sigma(v, i, Eps):
        E, nu = material_parameters[i]
        lmbda = E*nu/(1+nu)/(1-2*nu)
        mu = E/2/(1+nu)
        return lmbda*df.tr(eps(v) + Eps) * df.Identity(2) + 2*mu * (eps(v) + Eps)

    Ve = df.VectorElement("CG", mesh.ufl_cell(), 2)
    Re = df.VectorElement("R", mesh.ufl_cell(), 0)
    W = df.FunctionSpace(mesh, df.MixedElement([Ve, Re]), constrained_domain=PeriodicBoundary(mesh_corners, tolerance=1e-4))
    V = df.FunctionSpace(mesh, Ve)

    v_, lamb_ = df.TestFunctions(W)
    dv, dlamb = df.TrialFunctions(W)
    w = df.Function(W)
    N = df.FacetNormal(mesh)

    dx = df.Measure('dx', domain=mesh)
    # dS = df.Measure("ds", domain=mesh, subdomain_data=facet)

    Eps = df.Constant(((0, 0), (0, 0)))
    # internal pressure

    # Pressure is applied on physical surfaces 2 dS(2), and 3, dS(3)
    F = df.inner(sigma(dv, 0, Eps), eps(v_)) * dx \
    # - df.inner(-P*N, v_) * dS(2)

    a, L = df.lhs(F), df.rhs(F)
    a += df.dot(lamb_, dv) * dx + df.dot(dlamb, v_) * dx

    def Voigt2strain(s):
        return np.array([[s[0]   , s[2]/2.],
                         [s[2]/2., s[1]   ]])

    def get_macro_strain(i):
        """returns the macroscopic strain for the 3 elementary load cases"""
        Eps_Voigt = np.zeros(3)
        Eps_Voigt[i] = 1
        return Voigt2strain(Eps_Voigt)

    def stress2Voigt(s):
        return df.as_vector([s[0,0], s[1,1], s[0,1]])

    integrated_tissue_area = df.assemble(df.Constant(1.0) * df.dx(mesh))

    xdmf_file_per = df.XDMFFile("perturbations.xdmf")
    xdmf_file_sol = df.XDMFFile("solutions.xdmf")
    xdmf_file_sol_1 = df.XDMFFile("solutions_1.xdmf")
    xdmf_file_sol_2 = df.XDMFFile("solutions_2.xdmf")
    xdmf_file_sol_3 = df.XDMFFile("solutions_3.xdmf")
    u_tot = df.Function(V)
    u_tot.rename('u', 'u')
    u_tot_1 = df.Function(V)
    u_tot_1.rename('u', 'u')
    u_tot_2 = df.Function(V)
    u_tot_2.rename('u', 'u')
    u_tot_3 = df.Function(V)
    u_tot_3.rename('u', 'u')
    u_tensile_strain = df.Function(V)
    u_tensile_strain.rename('u', 'u')

    u_uc = df.Function(V)

    Chom = np.zeros((3, 3))
    for (j, case) in enumerate(["Exx", "Eyy", "Exy"]):
        #print("Solving {} case...".format(case))
        macro_strain = get_macro_strain(j)
        Eps.assign(df.Constant(macro_strain))
        df.solve(a == L, w, [], solver_parameters={"linear_solver": "cg"})
        (v, lamb) = df.split(w)
        # xdmf_file_per.write(w, float(j))
        Sigma = np.zeros((3,))
        for k in range(3):
            Sigma[k] = df.assemble(stress2Voigt(sigma(v, 0, Eps))[k] * dx)/vol
        Chom[j, :] = Sigma
    print("Chom =" +str(Chom))

    lmbda_hom = Chom[0, 1]
    mu_hom = Chom[2, 2]
    print("lmbda_bar:", lmbda_hom)
    print("mu_bar:", mu_hom)

    E_hom = mu_hom*(3*lmbda_hom + 2*mu_hom)/(lmbda_hom + mu_hom)
    nu_hom = lmbda_hom/(lmbda_hom + mu_hom)/2
    print("E_bar:", E_hom)
    print("nu_bar:", nu_hom)

    print(Chom[0, 0], lmbda_hom + 2*mu_hom)

    print("C11:" +str(E_hom*(1 - nu_hom)/(1+nu_hom)/(1-2*nu_hom)))

    Shom = np.linalg.inv(Chom)
    print("Shom =" +str(Shom))

    E_hom = 1./Shom[0,0]
    print("E_bar:", E_hom)

    E_hom = 1./Shom[1,1]
    print("E_bar:", E_hom)

    nu_hom = -Shom[0,1]/Shom[0,0]
    print("nu_bar:", nu_hom)

    nu_hom = -Shom[1,0]/Shom[1,1]
    print("nu_bar:", nu_hom)

    G_hom = 1./Shom[2,2]
    print("G_bar:", G_hom)

    G_hom = E_hom/2/(1+nu_hom)
    print("G_bar:", G_hom)

    h = 0.25
    l = 1 - 2 * math.sqrt(3) * h/3
    porous_area = 3*math.sqrt(3)*l**2
    tissue_area = (vol - porous_area)

    calculated_porosity = (vol-integrated_tissue_area)/vol
    Phi_s_0 = integrated_tissue_area/vol

    print("Porosity = " + str(calculated_porosity))


   
    return([lmbda_hom, mu_hom, calculated_porosity])
