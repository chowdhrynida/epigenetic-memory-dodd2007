import numpy as np
import matplotlib.pyplot as plt

# Define parameters
N = 60  # Number of nucleosomes
F_values = [0.4, 1.0, 1.4, 2.0]  # Feedback-to-noise ratios
Tmax = 5000 * N  # Total time steps (average attempted conversions per nucleosome)
Teq = 10  # Equilibrium time steps

# Define nucleosome states
A, U, M = 0, 1, 2

# Function to update nucleosome state
def step(state, F):
    # Calculate recruitment probability
    alpha = F / (F + 1)
    
    # Choose a random nucleosome
    n1 = np.random.randint(N)
    
    # Recruited conversion
    if np.random.rand() < alpha:
        n2 = np.random.randint(N)
        if state[n2] == M:
            if state[n1] == A:
                state[n1] = U
            elif state[n1] == U:
                state[n1] = M
        elif state[n2] == A:
            if state[n1] == M:
                state[n1] = U
            elif state[n1] == U:
                state[n1] = A
    # Noisy conversion
    else:
        r = np.random.rand()
        if state[n1] == A:
            if r < 1/3:
                state[n1] = U
        elif state[n1] == U:
            if r < 1/3:
                state[n1] = A
            elif r < 2/3:
                state[n1] = M
        elif state[n1] == M:
            if r < 1/3:
                state[n1] = U
    return state

# Create a figure with subplots
fig, axs = plt.subplots(len(F_values), 2, figsize=(12, 8))

# Initialize list to store G values
G_values = []

# Iterate over F values
for idx, F in enumerate(F_values):
    # Initialize nucleosome states
    state = np.random.choice([A, M, U], size=N)
    
    # Equilibrium steps
    for t in range(Teq):
        state = step(state, F)
    
    # Production run
    M_counts = []
    U_counts = []
    A_counts = []
    G_record = []
    time = []
    t = 0
    for i in range(Tmax):
        state = step(state, F)
        t += 1 / N  # increment time by 1/N (average attempted conversions per nucleosome)
        if i % N == 0:
            # Record M count every N steps
            M1 = np.sum(state == M)
            A_count = np.sum(state == A)
            U_count = np.sum(state == U)
            if M1 + A_count > 0:
                G_record.append(np.abs(M1 - A_count) / (M1 + A_count))
            M_counts.append(np.sum(state == M))
            A_counts.append(np.sum(state == A))
            U_counts.append(np.sum(state == U))
            time.append(t)
    
    # Calculate G
    G = np.mean(G_record)
    G_values.append(G)
    print(f"F = {F}, G = {G:.2f}")
    
    # Plot M vs time
    axs[idx, 0].plot(time, M_counts)
    axs[idx, 0].set_title(f"F = {F}")
    axs[idx, 0].set_xlabel("Time")
    axs[idx, 0].set_ylabel("M")
    
    # Plot P(M) vs M
    axs[idx, 1].hist(M_counts, bins=np.arange(0, N+1), density=True, align='left', rwidth=0.8)
    axs[idx, 1].set_title(f"F = {F}")
    axs[idx, 1].set_xlabel("M")
    axs[idx, 1].set_ylabel("Probability")

# Layout so plots do not overlap
plt.tight_layout()

# Show the plot
plt.show()

# Calculate alpha values
for F in F_values:
    alpha = F / (F + 1)
    print(f"F = {F}, alpha = {alpha:.2f}")