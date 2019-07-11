import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from pandas import DataFrame


class SolutionObject(DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def fast_plot(self, x_col, y_col, z_col=None, fig_size=(9, 6),
                  *args, **kwargs):
        """
        The arguments x_col, y_col and z_col must be an integer or an array. If
        the argument is an integer then the data for that column will be takan
        as the column number corresponding to this integer (for example, if
        x_col=1 then the data for the x label will be the second column of the
        solution).

        If an array is passed, the data for that lebel will be that array.

        The plot will be a 2D plot but default, unless the data for the z
        label is passed.

        """

        data = []
        for coord in (x_col, y_col, z_col):
            if isinstance(coord, int):
                data.append(self[coord])
            else:
                data.append(coord)

        # 2D plot
        if z_col is None:
            fig, ax = plt.subplots(figsize=fig_size)
            ax.plot(data[0], data[1], *args, **kwargs)

            return ax

        # 3D plot
        fig = plt.figure(frameon=0)
        proj_3D = fig.gca(projection='3d')
        proj_3D.plot(data[0], data[1], data[2])

        return proj_3D


class DiffEqSolver:

    def __init__(self, dif_equation, method='rk4'):
        """

        dif_equation must take by the arguments the time t and vector y
        and return a numpy like array with the same length as
        y.

        """

        self.dy = dif_equation
        self.method = method
        self._solution = None

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, m):
        valid_methods = ['rk4', 'rk2', 'newton']

        try:
            valid_method = m.lower() in valid_methods
        except AttributeError:
            raise AttributeError("Argument method must be a string")

        if valid_method is False:
            msg = "Unknown method. Valid methods are: 'rk4', 'rk2', 'newton'"
            raise TypeError(msg)

        self._method = m

    def next_step_newton(self, t, y):
        raise AttributeError('This method is not implemented!')

    def next_step_rk2(self, t, y):
        k1 = self.dy(t, y)
        k2 = self.dy(t + self.h, y + self.h * k1)

        return y + (self.h / 2.0) * (k1 + k2)

    def next_step_rk4(self, t, y):

        k1 = self.dy(t, y)
        k2 = self.dy(t + 0.5 * self.h, y + 0.5 * self.h * k1)
        k3 = self.dy(t + 0.5 * self.h, y + 0.5 * self.h * k2)
        k4 = self.dy(t + self.h, y + self.h * k3)

        return y + (self.h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    def get_trajectory(self, y0, next_step):
        '''
        Here it is performed the main loop to find the trajectory

        '''

        # Number of steps
        n_steps = len(self.time_steps)

        # Initial conditions
        sol_eqs = np.empty((n_steps, len(y0)))
        sol_eqs[0] = y0

        # Loop for solving the equations
        for i, t in enumerate(self.time_steps[1:], 1):
            sol_eqs[i] = next_step(t, sol_eqs[i-1])

        return sol_eqs

    def solve_eqs(self, initial_conditions, ti, tf, h, index=None, column_names=None,
                  numpy_array=False, save=False, file_name="sol_eqs.csv"):

        # Initial conditions
        self.initial_conditions = initial_conditions
        self.initial_time, self.final_time = ti, tf
        self.h = h

        # Choose the method
        if self.method.lower() == 'rk4':
            next_step = self.next_step_rk4
        elif self.method.lower() == 'rk2':
            next_step = self.next_step_rk2
        elif self.method.lower() == 'newton':
            next_step = self.next_step_newton

        # Time vector
        n_steps = int((tf - ti) / h)
        self.time_steps = np.linspace(ti, tf, n_steps)

        # Get the trajectory:
        sol_eqs = self.get_trajectory(initial_conditions, next_step)

        # Save the file if save is True:
        if save is True:
            np.savetxt(file_name, sol_eqs)

        # Create a dataframe
        self._solution = SolutionObject(
            sol_eqs, index=None, columns=column_names)

        # Return a dataframe unless numpy_array is not False
        if numpy_array:
            return sol_eqs

        return self._solution
