import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
BR = 5
scale = 3
mean = 0
decay = -1/3
decay2 = -1/3

standard_deviation1 = 1.19/2
standard_deviation2 = 1.19/2
standard_deviation3 = 1.19/4
x_values = np.arange(-5, 5, .01)
normaldist1 = norm(mean, standard_deviation1)
y_values1 =  normaldist1.pdf(x_values)
y_values1 = ((scale-1)*y_values1/max(y_values1))
normaldist2 = norm(mean, standard_deviation2)
y_values2 = normaldist2.pdf(x_values)
y_values2 = ((scale)*y_values2/max(y_values2))

yesarray = ((x_values>-1.19) & (x_values<1.19))

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, sharey = True)
plt.setp((ax, ax2), xticks=[-4.8, -3.6,-2.4,-1.2,0, 1.2, 2.4, 3.6, 4.8],
        yticks=np.linspace(BR, BR*(scale+1), 7))


y_meas = BR*3

ax.fill_between(x_values,BR,y_meas, yesarray, color = 'red', alpha = .45, linewidths = 1,label = "measured")
ax2.fill_between(x_values,BR,y_meas, yesarray, color = 'red', alpha = .45, linewidths = 1,label = "measured")
#for i in range(0,10):
#    y_meas -= ((BR*3)-BR)/10
#    ax.fill_between(x_values,BR,y_meas, yesarray, color = 'red', alpha = .15, linewidths = 1,label = "_nolabel_")
   # ax2.fill_between(x_values,BR,y_meas, yesarray, color = 'red', alpha = .15, linewidths = 1,label = #"_nolabel_")
#    print(y_meas)
ax.plot(x_values,(y_values1+1)*BR, color = "black",label = "modeled decay")
print(np.mean((y_values1[yesarray]+1)*BR))
alph = 1
for i in range(1,12):
    alph -= .075
    y_values1 = y_values1*np.exp(decay2)
    ax.plot(x_values,(y_values1+1)*BR, color = "black",  alpha = alph,  linewidth = 1,label = "_nolegend_")
    print(np.mean((y_values1[yesarray]+1)*BR))
ax.set_ylabel('migration rate (m/year)')
ax.set_title("Case 1")
ax.legend(loc = "upper left")


ax2.plot(x_values,(y_values2+1)*BR, color = "black", label = "modeled decay")
print(np.mean((y_values2[yesarray]+1)*BR))
alph = 1
for i in range(1,12):
    alph -=.075
    y_values2 = y_values2*np.exp(decay)
    ax2.plot(x_values,(y_values2+1)*BR, color = "black", alpha = alph,  linewidth = 1, label = "_nolegend_")
    print(np.mean((y_values2[yesarray]+1)*BR))
ax2.set_xlabel('distance from cutoff bend (length removed by cutoff)')
ax2.set_ylabel('migration rate (m/year)')
ax2.set_title("Case 2")

plt.savefig("nonlocaleffectshape.png", dpi = 1500)
plt.close()