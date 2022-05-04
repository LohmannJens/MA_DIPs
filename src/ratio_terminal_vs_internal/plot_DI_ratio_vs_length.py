import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression

xlsx_path = os.path.join("data", "SchwartzLowen_Fig3a_MdckP3.xlsx")
df = pd.read_excel(xlsx_path, skiprows=[9,10,11,12])

x = np.array(df["Length"]).reshape((-1, 1))
y = np.array(df["average"])
err = np.array(df["standard_deviation"])

fig, ax = plt.subplots(figsize = (9, 9))
ax.scatter(x, y, s=60, alpha=0.7, edgecolors="k")
ax.errorbar(x, y, yerr=err, fmt="o", capsize=5.0)

#model = LinearRegression().fit(x, y)
zero_model = LinearRegression(fit_intercept=False).fit(x, y)

model_x = np.linspace(0, 2341, num=100).reshape((-1, 1))

#ax.plot(x_model, model.predict(x_model), label=f"{model.score(x, y):.4f}")
ax.plot(model_x, zero_model.predict(model_x), label="linear model starting at 0")

exp_x = x.reshape((1, -1))[0]
exp_model_x = model_x.reshape((1, -1))[0]
exp_model = np.polyfit(exp_x, np.log(y), 1)
ax.plot(exp_model_x, np.exp(exp_model[1]) * np.exp(exp_model[0] * exp_model_x), label="exponential model")


ax.set_xlim(left=0)
plt.legend(loc="upper left")


plt.show()
