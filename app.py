import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
import ctypes
import os

mylib = ctypes.CDLL("D:\\IMPORTANT DATA\\AIU Computer Engineering\\CE Year 2\\Semester 2\\MAT315 Numerical Analysis\Project\\Numerical-Methods-In-C\\impl.dll")


# // Root-finding Methods Section:
# // - Bisection Method
# // - Newton-Raphson Method
# // - Secant Method
# // - Regula Falsi Method (False Position)
# // - Fixed-Point Iteration Method
# // - Halley's Method

def bisection(f, a, b, err, iter):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    return mylib.bisection(func, a, b, err, iter)

def newt_raph(f, df, x0, err, iter):
    CFUNCTYPE_f = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    CFUNCTYPE_df = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func_f = CFUNCTYPE_f(f)
    func_df = CFUNCTYPE_df(df)
    return mylib.newt_raph(func_f, func_df, x0, err, iter)

def secant(f, x0, x1, err, iter):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    return mylib.secant(func, x0, x1, err, iter)

def reg_falsi(f, a, b, err, iter):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    return mylib.reg_falsi(func, a, b, err, iter)

def fixed_pt_iter(g, x0, err, iter):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(g)
    return mylib.fixed_pt_iter(func, x0, err, iter)

def halleys(f, df, ddf, x0, err, iter):
    CFUNCTYPE_f = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    CFUNCTYPE_df = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    CFUNCTYPE_ddf = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func_f = CFUNCTYPE_f(f)
    func_df = CFUNCTYPE_df(df)
    func_ddf = CFUNCTYPE_ddf(ddf)
    return mylib.halleys(func_f, func_df, func_ddf, x0, err, iter)

# // Interpolation & Approximation Methods Section:
# // - Lagrange Interpolation
# // - Divided Differences Interpolation
# // - Forward Differences Interpolation
# // - Backward Differences Interpolation

def lagrange(x, y, n, x_interp):
    x_arr = (ctypes.c_double * n)(*x)
    y_arr = (ctypes.c_double * n)(*y)
    mylib.lagrange.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double]
    mylib.lagrange.restype = ctypes.c_double
    return mylib.lagrange(x_arr, y_arr, n, x_interp)

def divi_diffs(x, y, n):
    x_arr = (ctypes.c_double * n)(*x)
    y_arr = (ctypes.c_double * n)(*y)
    coefficients = (ctypes.c_double * n)()
    mylib.divi_diffs.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    mylib.divi_diffs.restype = ctypes.c_double
    mylib.divi_diffs(x_arr, y_arr, coefficients, n)
    return list(coefficients)

def fwd_intpol(x, y, n, x_interp):
    x_arr = (ctypes.c_double * n)(*x)
    y_arr = (ctypes.c_double * n)(*y)
    mylib.fwd_intpol.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double]
    mylib.fwd_intpol.restype = ctypes.c_double
    return mylib.fwd_intpol(x_arr, y_arr, n, x_interp)

def bckwd_intpol(x, y, n, x_interp):
    x_arr = (ctypes.c_double * n)(*x)
    y_arr = (ctypes.c_double * n)(*y)
    mylib.bckwd_intpol.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double]
    mylib.bckwd_intpol.restype = ctypes.c_double
    return mylib.bckwd_intpol(x_arr, y_arr, n, x_interp)

# // Linear Algebra Methods Section:
# // - Jacobi Iterative Method
# // - Gauss-Seidel Iterative Method

def jacobi(A, b, x, n, err, iter):
    A_arr = (ctypes.POINTER(ctypes.c_double) * n)()
    for i in range(n):
        A_arr[i] = (ctypes.c_double * n)(*A[i])
    b_arr = (ctypes.c_double * n)(*b)
    x_arr = (ctypes.c_double * n)(*x)
    mylib.jacobi.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.c_int]
    mylib.jacobi.restype = None
    mylib.jacobi(A_arr, b_arr, x_arr, n, err, iter)
    return list(x_arr)

def gau_sei(A, b, x, n, err, iter):
    A_arr = (ctypes.POINTER(ctypes.c_double) * n)()
    for i in range(n):
        A_arr[i] = (ctypes.c_double * n)(*A[i])
    b_arr = (ctypes.c_double * n)(*b)
    x_arr = (ctypes.c_double * n)(*x)
    mylib.gau_sei.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.c_int]
    mylib.gau_sei.restype = None
    mylib.gau_sei(A_arr, b_arr, x_arr, n, err, iter)
    return list(x_arr)

# // Numerical Differentiation Methods Section:
# // - Two Point Forward Difference Method
# // - Two Point Backward Difference Method
# // - Three Point Forward Difference Method
# // - Three Point Backward Difference Method
# // - Three Point Central Difference Method

def two_pt_fwd_diff(f, x, h):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    mylib.two_pt_fwd_diff.argtypes = [CFUNCTYPE, ctypes.c_double, ctypes.c_double]
    mylib.two_pt_fwd_diff.restype = ctypes.c_double
    return mylib.two_pt_fwd_diff(func, x, h)

def two_pt_bckd_diff(f, x, h):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    mylib.two_pt_bckd_diff.argtypes = [CFUNCTYPE, ctypes.c_double, ctypes.c_double]
    mylib.two_pt_bckd_diff.restype = ctypes.c_double
    return mylib.two_pt_bckd_diff(func, x, h)

def three_pt_fwd_diff(f, x, h):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    mylib.three_pt_fwd_diff.argtypes = [CFUNCTYPE, ctypes.c_double, ctypes.c_double]
    mylib.three_pt_fwd_diff.restype = ctypes.c_double
    return mylib.three_pt_fwd_diff(func, x, h)

def three_pt_bckd_diff(f, x, h):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    mylib.three_pt_bckd_diff.argtypes = [CFUNCTYPE, ctypes.c_double, ctypes.c_double]
    mylib.three_pt_bckd_diff.restype = ctypes.c_double
    return mylib.three_pt_bckd_diff(func, x, h)

def three_pt_cent_diff(f, x, h):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    mylib.three_pt_cent_diff.argtypes = [CFUNCTYPE, ctypes.c_double, ctypes.c_double]
    mylib.three_pt_cent_diff.restype = ctypes.c_double
    return mylib.three_pt_cent_diff(func, x, h)

# // Numerical Integration Methods Section:
# // - Trapezoidal Rule
# // - Simpson's Rule

def trap(f, a, b, n):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    mylib.trap.argtypes = [CFUNCTYPE, ctypes.c_double, ctypes.c_double, ctypes.c_int]
    mylib.trap.restype = ctypes.c_double
    return mylib.trap(func, a, b, n)

def simps(f, a, b, n):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    mylib.simps.argtypes = [CFUNCTYPE, ctypes.c_double, ctypes.c_double, ctypes.c_int]
    mylib.simps.restype = ctypes.c_double
    return mylib.simps(func, a, b, n)

# // Ordinary Differential Equations Section:
# // - Euler's Method
# // - Modified Euler's Method
# // - Runge-Kutta's 2nd Order
# // - Runge-Kutta's 3rd Order
# // - Runge-Kutta's 4th Order

def euler(f, x0, y0, h, n):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    x = (ctypes.c_double * (n + 1))()
    y = (ctypes.c_double * (n + 1))()
    mylib.euler(func, x0, y0, h, n, x, y)
    return list(x), list(y)

def mod_euler(f, x0, y0, h, n):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    x = (ctypes.c_double * (n + 1))()
    y = (ctypes.c_double * (n + 1))()
    mylib.mod_euler(func, x0, y0, h, n, x, y)
    return list(x), list(y)

def rk2(f, x0, y0, h, n):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    x = (ctypes.c_double * (n + 1))()
    y = (ctypes.c_double * (n + 1))()
    mylib.rk2(func, x0, y0, h, n, x, y)
    return list(x), list(y)

def rk3(f, x0, y0, h, n):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    x = (ctypes.c_double * (n + 1))()
    y = (ctypes.c_double * (n + 1))()
    mylib.rk3(func, x0, y0, h, n, x, y)
    return list(x), list(y)

def rk4(f, x0, y0, h, n):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    x = (ctypes.c_double * (n + 1))()
    y = (ctypes.c_double * (n + 1))()
    mylib.rk4(func, x0, y0, h, n, x, y)
    return list(x), list(y)


# // Aitken's Delta-Squared Acceleration Process
def atkn_dlta_sqd(x):
    n = len(x)
    array_type = ctypes.c_double * n
    mylib.atkn_dlta_sqd.argtypes = [array_type, ctypes.c_int]
    mylib.atkn_dlta_sqd.restype = ctypes.c_double
    
    x_array = array_type(*x)
    return mylib.atkn_dlta_sqd(x_array, n)

class NumericalMethodsApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Numerical Methods")
        self.geometry(f"{self.winfo_screenwidth()}x{self.winfo_screenheight()}")
        
        # Set the icon for the window
        if os.name == 'nt':  # Windows
            self.iconbitmap('D:\\IMPORTANT DATA\\AIU Computer Engineering\\CE Year 2\\Semester 2\\MAT315 Numerical Analysis\\Project\\Numerical-Methods-In-C\\nummethodslogoJPG.ico')

        # Create sections
        self.create_root_finding_section()
        self.create_interpolation_section()
        self.create_linear_alg_section()
        self.create_num_diff_section()
        self.create_num_integ_section()
        self.create_ode_section()

    def create_root_finding_section(self):
        frame = tk.Frame(self)
        frame.pack(pady=10)

        label = tk.Label(frame, text="Root-Finding Methods", font=("Helvetica", 16))
        label.pack()

        bisection_button = tk.Button(frame, text="Bisection Method", command=self.bisection_method)
        bisection_button.pack(pady=5)
        
        # Add more buttons for other root-finding methods...

    def create_interpolation_section(self):
        frame = tk.Frame(self)
        frame.pack(pady=10)

        label = tk.Label(frame, text="Interpolation Methods", font=("Helvetica", 16))
        label.pack()

        # Add buttons for interpolation methods...

    def create_linear_alg_section(self):

        frame = tk.Frame(self)
        frame.pack(pady=10)

        label = tk.Label(frame, text="Linear Algebra Methods", font=("Helvetica", 16))
        label.pack()

        # Add buttons for linear algebra methods...

    def create_num_diff_section(self):
        frame = tk.Frame(self)
        frame.pack(pady=10)

        label = tk.Label(frame, text="Numerical Differentiation Methods", font=("Helvetica", 16))
        label.pack()

        # Add buttons for numerical differentiation methods...

    def create_num_integ_section(self):
        frame = tk.Frame(self)
        frame.pack(pady=10)

        label = tk.Label(frame, text="Numerical Integration Methods", font=("Helvetica", 16))
        label.pack()

        # Add buttons for numerical integration methods...

    def create_ode_section(self):
        frame = tk.Frame(self)
        frame.pack(pady=10)

        label = tk.Label(frame, text="Ordinary Differential Equation Methods", font=("Helvetica", 16))
        label.pack()

        # Add buttons for ODE methods...

    # Define command methods for buttons
    def bisection_method(self):
        # Example popup to get inputs
        popup = tk.Toplevel(self)
        popup.title("Bisection Method")
        
        tk.Label(popup, text="Function:").grid(row=0, column=0)
        function_entry = tk.Entry(popup)
        function_entry.grid(row=0, column=1)

        tk.Label(popup, text="a:").grid(row=1, column=0)
        a_entry = tk.Entry(popup)
        a_entry.grid(row=1, column=1)
        
        tk.Label(popup, text="b:").grid(row=2, column=0)
        b_entry = tk.Entry(popup)
        b_entry.grid(row=2, column=1)
        
        tk.Label(popup, text="Tolerance:").grid(row=3, column=0)
        tol_entry = tk.Entry(popup)
        tol_entry.grid(row=3, column=1)
        
        tk.Label(popup, text="Max Iterations:").grid(row=4, column=0)
        max_iter_entry = tk.Entry(popup)
        max_iter_entry.grid(row=4, column=1)
        
        def run_bisection():
            try:
                f = lambda x: eval(function_entry.get())
                a = float(a_entry.get())
                b = float(b_entry.get())
                tol = float(tol_entry.get())
                max_iter = int(max_iter_entry.get())
                
                result = bisection(f, a, b, tol, max_iter)
                messagebox.showinfo("Result", f"Root: {result}")
            except Exception as e:
                messagebox.showerror("Error", str(e))

        tk.Button(popup, text="Run", command=run_bisection).grid(row=5, columnspan=2)


# Run the app
if __name__ == "__main__":
    app = NumericalMethodsApp()
    app.mainloop()

# class NumericalMethodsApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Numerical Methods")
        self.geometry(f"{self.winfo_screenwidth()}x{self.winfo_screenheight()}")

        # Set the icon for the window
        if os.name == 'nt':
            self.iconbitmap('D:\\IMPORTANT DATA\\AIU Computer Engineering\\CE Year 2\\Semester 2\\MAT315 Numerical Analysis\\Project\\Numerical-Methods-In-C\\nummethodslogoJPG.ico')

        self.frames = {}
        container = ttk.Frame(self)
        container.pack(fill="both", expand=True)

        for F in (HomePage, RootFindingPage, InterpolationPage, LinearAlgebraPage,
                  NumericalIntegrationPage, NumericalDifferentiationPage, ODEsPage):
            page_name = F.__name__
            frame = F(parent=container, controller=self)
            self.frames[page_name] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame("HomePage")

    def show_frame(self, page_name):
        frame = self.frames[page_name]
        frame.tkraise()

# class HomePage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Numerical Methods", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        buttons = [
            ("Root-finding methods", "RootFindingPage"),
            ("Interpolation and approximation", "InterpolationPage"),
            ("Linear Algebra", "LinearAlgebraPage"),
            ("Numerical Integration", "NumericalIntegrationPage"),
            ("Numerical Differentiation", "NumericalDifferentiationPage"),
            ("ODEs", "ODEsPage")
        ]

        for text, page_name in buttons:
            button = ttk.Button(self, text=text, command=lambda pn=page_name: controller.show_frame(pn))
            button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

# class RootFindingPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Root-finding Methods", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        buttons = [
            ("Bisection Method", "BisectionPage"),
            ("Newton-Raphson Method", "NewtonRaphsonPage"),
            ("Secant Method", "SecantPage"),
            ("False Position Method", "FalsePositionPage"),
            ("Fixed-Point Iteration Method", "FixedPointIterationPage"),
            ("Halley's Method", "HalleysMethodPage")
        ]

        for text, page_name in buttons:
            button = ttk.Button(self, text=text, command=lambda pn=page_name: controller.show_frame(pn))
            button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

        back_button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        back_button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

# class InterpolationPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Interpolation and Approximation", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        buttons = [
            ("Lagrange Interpolation", "LagrangeInterpolationPage"),
            ("Divided Differences", "DividedDifferencesPage"),
            ("Forward Interpolation", "ForwardInterpolationPage"),
            ("Backward Interpolation", "BackwardInterpolationPage")
        ]

        for text, page_name in buttons:
            button = ttk.Button(self, text=text, command=lambda pn=page_name: controller.show_frame(pn))
            button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

        back_button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        back_button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

# class LinearAlgebraPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Linear Algebra", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        buttons = [
            ("Jacobi Method", "JacobiMethodPage"),
            ("Gauss-Seidel Method", "GaussSeidelMethodPage")
        ]

        for text, page_name in buttons:
            button = ttk.Button(self, text=text, command=lambda pn=page_name: controller.show_frame(pn))
            button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

        back_button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        back_button.pack(side="top", pady=5, anchor="center")

# class NumericalIntegrationPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Numerical Integration", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        buttons = [
            ("Trapezoidal Rule", "TrapezoidalRulePage"),
            ("Simpson's Rule", "SimpsonsRulePage")
        ]

        for text, page_name in buttons:
            button = ttk.Button(self, text=text, command=lambda pn=page_name: controller.show_frame(pn))
            button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

        back_button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        back_button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

# class NumericalDifferentiationPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Numerical Differentiation", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        buttons = [
            ("2 Points Formula (Forward)", "TwoPointsForwardFormulaPage"),
            ("2 Points Formula (Backward)", "BackwardTwoPointsFormulaPage"),
            ("3 Points Formula (Forward)", "ThreePointsForwardFormulaPage"),
            ("3 Points Formula (Backward)", "ThreePointsBackwardFormulaPage"),
            ("3 Points Formula (Central)", "ThreePointsCentralFormulaPage")
        ]

        for text, page_name in buttons:
            button = ttk.Button(self, text=text, command=lambda pn=page_name: controller.show_frame(pn))
            button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

        back_button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        back_button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

# class ODEsPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Ordinary Differential Equations (ODEs)", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        buttons = [
            ("Euler's Method", "EulersMethodPage"),
            ("Modified Euler's Method", "ModifiedEulersMethodPage"),
            ("Runge-Kutta (RK2)", "RungeKutta2Page"),
            ("Runge-Kutta (RK3)", "RungeKutta3Page"),
            ("Runge-Kutta (RK4)", "RungeKutta4Page")
        ]

        for text, page_name in buttons:
            button = ttk.Button(self, text=text, command=lambda pn=page_name: controller.show_frame(pn))
            button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")

        back_button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        back_button.pack(side="top", pady=5, anchor="center", expand=False, fill="none")
















































































































































