import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
import ctypes
import ctypes
import os

# if os.name == 'nt':
#     mylib = ctypes.CDLL("build/impl.dll")

def bisection(f, a, b, err, iter):
    CFUNCTYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
    func = CFUNCTYPE(f)
    return mylib.bisection(func, a, b, err, iter)

# class NumericalMethodsApp(tk.Tk):
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

class NumericalMethodsApp(tk.Tk):
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

class HomePage(ttk.Frame):
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

class RootFindingPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Root-finding Methods", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        button.pack(side="top", pady=5, anchor="center")


class InterpolationPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Interpolation and Approximation", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        button.pack(side="top", pady=5, anchor="center")

class LinearAlgebraPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Linear Algebra", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        button.pack(side="top", pady=5, anchor="center")

class NumericalIntegrationPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Numerical Integration", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        button.pack(side="top", pady=5, anchor="center")

class NumericalDifferentiationPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Numerical Differentiation", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        button.pack(side="top", pady=5, anchor="center")

class ODEsPage(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        label = ttk.Label(self, text="Ordinary Differential Equations (ODEs)", font=("Helvetica", 16))
        label.pack(pady=10, padx=10)

        button = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame("HomePage"))
        button.pack(side="top", pady=5, anchor="center")
# Run the app
if __name__ == "__main__":
    app = NumericalMethodsApp()
    app.mainloop()