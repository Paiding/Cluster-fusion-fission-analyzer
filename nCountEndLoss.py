import numpy as np

sum_routes = np.zeros(1000)
for i in range(0,5000):
    with open(f'log/conf{i}_gain.dat', 'r') as file:
        lines = file.readlines()
        numbers = [float(line.strip()) for line in lines]
        for j in range(0,1000):
            sum_routes[j] = sum_routes[j] + numbers[j]

with open('sumgain.txt', 'w') as avg_file:
    for avg in sum_routes:
        avg_file.write(f'{avg}\n')
