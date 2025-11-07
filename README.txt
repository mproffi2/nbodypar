# Parallel N-Body Simulation

1. Compile:
   make nbody

2. Run Solar System simulation:
   ./nbody planet 200 5000000 10000 8 > solar.tsv

3. Run random 100 particles:
   ./nbody 100 1 10000 100 8 > random100.tsv

4. Visualize (Python script provided):
   python3 plot.py solar.tsv solar.pdf 1000

5. Try different thread counts (1, 2, 4, 8) and compare timings.