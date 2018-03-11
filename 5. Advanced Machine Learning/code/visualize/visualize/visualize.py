import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('HMMclass0.csv')
# df = pd.read_csv('HMMclass0.csv')
plt.plot(df['Klist'],df['train'],label='Train set')
plt.plot(df['Klist'],df['test'],label='Validation set')
plt.title('EM LOO cross-validation: Famous face')
plt.grid(linestyle=':')
plt.xlabel('States',fontsize=12)
plt.ylabel('loglike',fontsize=12)
plt.legend()
plt.show()
