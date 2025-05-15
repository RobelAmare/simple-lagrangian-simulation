import math  # Import the math module for mathematical functions
import matplotlib.pyplot as plt  # Import matplotlib for plotting
# Removed matplotlib.animation and numpy as they are no longer needed for this feature


def calculate_small_angle_period(l, g):  # Function to calculate the period using small angle approximation
    """
    Calculates the period of a simple pendulum using the small angle approximation.

    Args:
        l (float): Length of the pendulum in meters.
        g (float): Acceleration due to gravity in m/s^2.

    Returns:
        float: The period of the pendulum in seconds, or None if inputs are invalid.
    """
    if l <= 0 or g <= 0:  # Check for invalid (non-positive) input values
        # In a non-interactive environment, print errors and return None
        print("Error: Length and gravity must be positive values.")
        return None
    # Formula for the period of a simple pendulum (small angle approximation): T = 2 * pi * sqrt(l / g)
    period = 2 * math.pi * math.sqrt(l / g)  # Calculate the period
    return period  # Return the calculated period


def simulate_pendulum_motion(m, l, g, initial_theta, initial_theta_dot, duration, dt):  # Function to simulate pendulum motion
    """
    Simulates the motion of a simple pendulum using numerical integration (Euler-Cromer method).
    Also calculates kinetic, potential, and total energy at each step.

    Args:
        m (float): Mass of the pendulum bob in kg.
        l (float): Length of the pendulum in meters.
        g (float): Acceleration due to gravity in m/s^2.
        initial_theta (float): Initial angular position in radians.
        initial_theta_dot (float): Initial angular velocity in radians/s.
        duration (float): Total simulation time in seconds.
        dt (float): Time step for simulation in seconds.

    Returns:
        list: A list of tuples, where each tuple contains (time, angle, angular_velocity, kinetic_energy, potential_energy, total_energy).
              Returns an empty list if inputs are invalid.
    """
    if m <= 0 or l <= 0 or g <= 0 or duration <= 0 or dt <= 0:  # Check for invalid (non-positive) input values
         # Keep error message for invalid inputs
         print("Error: Mass, length, gravity, duration, and time step must be positive values.")
         return []
    if dt > duration:  # Check if time step is greater than duration
        # Keep error message for invalid inputs
        print("Error: Time step cannot be greater than simulation duration.")
        return []

    # Initialize variables
    time = 0.0  # Start time at 0
    theta = initial_theta  # Set initial angle
    theta_dot = initial_theta_dot  # Set initial angular velocity

    results = []  # List to store simulation results

    # Store initial state and energy
    kinetic_energy = 0.5 * m * (l * theta_dot)**2  # Calculate initial kinetic energy
    potential_energy = -m * g * l * math.cos(theta)  # Calculate initial potential energy
    total_energy = kinetic_energy + potential_energy  # Calculate initial total energy
    results.append((time, theta, theta_dot, kinetic_energy, potential_energy, total_energy))  # Store initial state


    # Numerical integration using Euler-Cromer method
    # The equation of motion is: theta_ddot = -(g/l) * sin(theta)
    # Use a small tolerance to ensure the last step is included if time is very close to duration
    while time <= duration + 1e-9:  # Loop until the end of the simulation duration
        # Calculate angular acceleration
        theta_ddot = -(g / l) * math.sin(theta)  # Compute angular acceleration

        # Update angular velocity using the acceleration
        new_theta_dot = theta_dot + theta_ddot * dt  # Update angular velocity

        # Update angle using the new angular velocity (Euler-Cromer)
        new_theta = theta + new_theta_dot * dt  # Update angle

        theta_dot = new_theta_dot  # Set new angular velocity
        theta = new_theta  # Set new angle
        time += dt  # Increment time

        # Calculate energy at the new state
        kinetic_energy = 0.5 * m * (l * theta_dot)**2  # Calculate kinetic energy
        potential_energy = -m * g * l * math.cos(theta)  # Calculate potential energy
        total_energy = kinetic_energy + potential_energy  # Calculate total energy

        results.append((time, theta, theta_dot, kinetic_energy, potential_energy, total_energy))  # Store results

    return results  # Return the simulation results


def plot_energy(simulation_results):  # Function to plot energy over time
    """
    Plots the Kinetic, Potential, and Total Energy over time using matplotlib.

    Args:
        simulation_results (list): A list of tuples from simulate_pendulum_motion.
    """
    if not simulation_results:  # Check if simulation results are empty
        print("No simulation data to plot energy.")
        return

    # Extract data for plotting
    times = [result[0] for result in simulation_results]  # Extract time values
    kinetic_energies = [result[3] for result in simulation_results]  # Extract kinetic energy values
    potential_energies = [result[4] for result in simulation_results]  # Extract potential energy values
    total_energies = [result[5] for result in simulation_results]  # Extract total energy values

    # Create the plot
    plt.figure(figsize=(10, 6)) # Set figure size for better readability
    plt.plot(times, kinetic_energies, label='Kinetic Energy', color='blue')  # Plot kinetic energy
    plt.plot(times, potential_energies, label='Potential Energy', color='red')  # Plot potential energy
    plt.plot(times, total_energies, label='Total Energy', color='green', linestyle='--') # Use dashed line for total energy

    # Add labels, title, and legend
    plt.xlabel("Time (s)")  # X-axis label
    plt.ylabel("Energy (J)")  # Y-axis label
    plt.title("Simple Pendulum Energy over Time")  # Plot title
    plt.legend()  # Show legend
    plt.grid(True)  # Show grid

    # Show the plot
    plt.show()  # Display the plot


# --- Command Line Interface ---

if __name__ == "__main__":  # Main entry point for script execution
    # Note: Interactive input using input() may not work in all execution environments.
    # For best results, run this script in a local terminal or command prompt.
    print("--- Simple Pendulum Numerical Analysis (User Input) ---")  # Print script title
    print("Please enter the requested values when prompted.")  # Prompt user for input
    print("Note: Interactive input may not work in all environments. Run locally for best results.")

    # Get user inputs
    try:
        m = float(input("Enter the mass of the pendulum bob (kg): "))  # Get mass from user
        l = float(input("Enter the length of the pendulum (m): "))  # Get length from user
        g = float(input("Enter the acceleration due to gravity (m/s^2): "))  # Get gravity from user
        initial_theta_deg = float(input("Enter the initial angular position (degrees): "))  # Get initial angle in degrees
        initial_theta_dot = float(input("Enter the initial angular velocity (rad/s): "))  # Get initial angular velocity
        duration = float(input("Enter the simulation duration (seconds): "))  # Get simulation duration
        dt = float(input("Enter the time step for simulation (seconds): "))  # Get time step

        # Convert initial angle from degrees to radians
        initial_theta = math.radians(initial_theta_deg)  # Convert degrees to radians

    except ValueError:
        print("\nInvalid input. Please enter numerical values.")  # Handle invalid input
        # Exit the script if input is invalid
        exit()
    except EOFError:
        # Handle EOFError specifically for environments that don't support input()
        print("\nError: Interactive input is not supported in this environment.")
        print("Please run this script in a local terminal or command prompt.")
        # Exit the script if input is not supported
        exit()


    print("\n--- Analysis Results ---")  # Print analysis results header

    # Calculate and display the period (small angle approximation)
    period = calculate_small_angle_period(l, g)  # Calculate period
    if period is not None:
        print(f"Calculated Period (small angle approximation): {period:.4f} seconds")  # Print period

    print("\n--- Simulation Results ---")  # Print simulation results header

    # Simulate motion and get results
    simulation_results = simulate_pendulum_motion(m, l, g, initial_theta, initial_theta_dot, duration, dt)  # Run simulation

    if simulation_results:
        # Print header for simulation results table
        # Adjusted column widths slightly for better alignment with potential larger numbers
        print(f"{'Time (s)':<10} | {'Angle (rad)':<15} | {'Angular Velocity (rad/s)':<22} | {'Kinetic Energy (J)':<18} | {'Potential Energy (J)':<20} | {'Total Energy (J)':<16}")
        print("-" * 10 + "-|-" + "-" * 15 + "-|-" + "-" * 22 + "-|-" + "-" * 18 + "-|-" + "-" * 20 + "-|-" + "-" * 16)


        # Print simulation data
        # Limit the number of output rows for very long simulations
        max_output_rows = 100  # Set maximum number of output rows
        output_interval = max(1, len(simulation_results) // max_output_rows)  # Calculate output interval

        for i, result in enumerate(simulation_results):  # Loop through simulation results
            if i % output_interval == 0:
                time, theta, theta_dot, ke, pe, te = result  # Unpack result tuple
                print(f"{time:<10.2f} | {theta:<15.4f} | {theta_dot:<22.4f} | {ke:<18.4f} | {pe:<20.4f} | {te:<16.4f}")  # Print result row

        if len(simulation_results) > max_output_rows:
            print(f"\n... Output truncated. Displaying {max_output_rows} rows out of {len(simulation_results)}. Consider a shorter duration or larger time step.")  # Notify if output is truncated

    print("\n----------------------------------------")  # Print separator

    # --- Energy Conservation Description ---
    print("\n--- Energy Conservation ---")  # Print energy conservation header
    print("In an ideal simple pendulum system (without friction or air resistance),")  # Print description
    print("mechanical energy (the sum of kinetic and potential energy) is conserved.")
    print("As the pendulum swings, energy is continuously converted between kinetic (energy of motion)")
    print("and potential (stored energy due to position in a gravitational field).")
    print("The total energy should remain constant over time, although numerical methods")
    print("may show slight variations due to approximation errors.")
    print("----------------------------------------")  # Print separator


    # --- Optional Energy Plotting ---
    # Ask the user if they want to see the energy plot
    if simulation_results: # Only ask if simulation was successful
        show_plot = input("Show energy plot? (yes/no): ").lower()  # Ask user if they want to see the plot
        if show_plot == 'yes':
            try:
                plot_energy(simulation_results)  # Plot energy if user agrees
            except Exception as e:
                print(f"\nError generating plot: {e}")  # Handle plotting errors
                print("Please ensure you have matplotlib installed (`pip install matplotlib`) and a compatible backend.")

