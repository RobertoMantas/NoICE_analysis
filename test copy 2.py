import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


# Some example data to display
x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)


fig, (ax1,ax2,ax3,ax4) = plt.subplots(4, sharex=True)
fig.suptitle('Vertically stacked subplots')
ax1.plot(x, y, label='Wind_Speed WRF')
ax2.plot(x, -y, label='Wind_Speed WRF')
ax3.plot(x, -y, label='Wind_Speed WRF')
ax4.plot(x, -y, label='Wind_Speed WRF')
ax4.set_ylabel('Wind Speed (m/s^2)',weight='bold')
ax4.set_xlabel('Time (m/s^2)',weight='bold')
ax1.legend(loc = 'upper left')
ax2.legend(loc = 'upper left')
ax3.legend(loc = 'upper left')
ax4.legend(loc = 'upper left')

plt.show()