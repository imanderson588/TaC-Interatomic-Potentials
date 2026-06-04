import matplotlib.pyplot as plt
import numpy as np
import yaml

with open("energy.yaml") as f:
    data = yaml.load(f, Loader=yaml.FullLoader)
fig, ax = plt.subplots()
for i in range(0,len(data)):
    ax.scatter(data[i]["stress"], data[i]["energy"],c='black')
plt.ylabel("Total Energy")
plt.xlabel("Stress")
plt.title("Peierls stress TaC Edge Dislocation")
plt.savefig("data-plot.png", dpi=300)
plt.show()
