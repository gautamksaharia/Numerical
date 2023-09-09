def run_simation(stochastic_process, N_sim, N_time):
    """
    simulate a stochastic process and return an array of results

    Parameters
    ----------
    stochastic process : ornstein process or random walk
        .
    N_sim : TYPE : int
        No of simualtion.
    N_time : float
        No of time points.

    Returns : an array of size (N_sim, N_time) containing the results of simulation
    -------
    None.

    """

  # initialize zero matrix of size (N_sim, N_time)
    all_simulation = np.zeros((N_sim, N_time))
    
    for i in range(N_sim):
        for j in range(1, N_time):
            all_simulation[i,j] = process(all_simulation[i, j-1])
    return all_simulation


ornstein_uhlenback= lambda x: x - 0.005*x + random.gauss(0, 1)
out = run_simation(ornstein_uhlenback, N_sim=100, N_time=100)

# Plot of random process
plt.plot(out.T)
plt.xlabel("time")
plt.ylabel("value")
plt.show()
