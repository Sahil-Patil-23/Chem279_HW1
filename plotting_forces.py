import numpy as np
import matplotlib.pyplot as plt 

def calculate_error(analytical_forces, numerical_forces):
    # Function to calculate the error that will be plotted
    error = 0.0
    for i in range(len(numerical_forces)):
        for j in range(3): # For each coordinate in 3D space (x,y,z)
            error += (numerical_forces[i][j] - analytical_forces[i][j]) ** 2
    return np.sqrt(error)


# Arrays below are all the data generated from 'HW1.cpp' file
analytical_forces = np.array([
    [-4.2191e-02 , -8.9489e-01 , -4.2191e-02, 9.7927e-01],
    [5.6398e-01 , 0.0000e+00 , -5.6398e-01 , 0.0000e+00],
    [0.0000e+00 , 0.0000e+00 , 0.0000e+00 , 0.0000e+00]
])

forward_forces_h_1 = np.array([
    [-4.6152e-02 , -8.4460e-01 , -4.6152e-02 , 1.0489e+00],
    [6.0210e-01 , 5.9650e-02 , -5.2940e-01 , -6.6810e-03],
    [-5.6147e-03 , -1.9989e-02 , -5.6147e-03 , -1.0635e-02]
])

central_forces_h_1 = np.array([
    [-4.2210e-02, -8.9817e-01, -4.2210e-02, 9.8259e-01],
    [5.6575e-01, -0.0000e+00, -5.6575e-01, 1.1102e-15],
    [-0.0000e+00, -0.0000e+00, -0.0000e+00, -0.0000e+00]
])

forward_forces_h_01 = np.array([
    [-4.2586e-02, -8.8957e-01, -4.2586e-02, 9.8592e-01],
    [5.6763e-01, 5.9503e-03, -5.6037e-01,-6.6903e-04],
    [-5.6188e-04, -2.0006e-03, -5.6188e-04, -1.0644e-03]
])

central_forces_h_01 = np.array([
    [-4.2191e-02, -8.9492e-01, -4.2191e-02, 9.7930e-01],
    [5.6400e-01, -0.0000e+00, -5.6400e-01, -1.1102e-14],
    [-0.0000e+00, -0.0000e+00, -0.0000e+00, -0.0000e+00]
])

forward_forces_h_001 = np.array([
    [-4.2231e-02, -8.9435e-01, -4.2231e-02, 9.7993e-01],
    [5.6435e-01, 5.9502e-04, -5.6362e-01, -6.6904e-05],
    [-5.6189e-05, -2.0006e-04, -5.6189e-05, -1.0644e-04]
])

central_forces_h_001 = np.array([
    [-4.2191e-02, -8.9489e-01, -4.2191e-02, 9.7927e-01],
    [5.6398e-01, -0.0000e+00, -5.6398e-01, -0.0000e+00],
    [-0.0000e+00, -0.0000e+00, -0.0000e+00, -0.0000e+00]
])

forward_forces_h_0001 = np.array([
    [-4.2195e-02, -8.9483e-01, -4.2195e-02, 9.7933e-01],
    [5.6402e-01, 5.9502e-05, -5.6395e-01, -6.6904e-06],
    [-5.6189e-06, -2.0006e-05, -5.6189e-06, -1.0644e-05]
])

central_forces_h_0001 = np.array([
    [-4.2191e-02, -8.9489e-01, -4.2191e-02, 9.7927e-01],
    [5.6398e-01, -0.0000e+00, -5.6398e-01, -1.1102e-12],
    [-0.0000e+00, -0.0000e+00, -0.0000e+00, -0.0000e+00]
])



h = np.array([0.1, 0.01, 0.001, 0.0001])

# Calculating error for forward differences
forward_errors = np.array([
    calculate_error(analytical_forces, forward_forces_h_1),
    calculate_error(analytical_forces, forward_forces_h_01),
    calculate_error(analytical_forces, forward_forces_h_001),
    calculate_error(analytical_forces, forward_forces_h_0001)
])

# Calculating error for central differences
central_errors = np.array([
    calculate_error(analytical_forces, central_forces_h_1),
    calculate_error(analytical_forces, central_forces_h_01),
    calculate_error(analytical_forces, central_forces_h_001),
    calculate_error(analytical_forces, central_forces_h_0001)
])

# Replace zeros in central_errors with a small value (e.g., 1e-10)
central_errors[central_errors == 0] = 1e-10

# Create a log-log plot
plt.figure(figsize=(10, 6))

slope_forward_d = np.polyfit(np.log(h), np.log(forward_errors), 1)[0]
slope_central_d = np.polyfit(np.log(h), np.log(central_errors), 1)[0]

# Forward difference error
plt.plot(np.log(h), np.log(forward_errors), label=f"Forward Difference Error (slope = {slope_forward_d})", marker='o')

# Central difference error
plt.plot(np.log(h), np.log(central_errors), label=f"Central Difference Error (slope = {slope_central_d})", marker='o')

# Labels and legend
plt.xlabel("Log(Step size h)")
plt.ylabel("Log(Error)")
plt.title("Error vs Step Size")
plt.legend()
plt.grid(True)

# Save the plot as a png file 
plt.savefig('Error_vs_StepSize.png')