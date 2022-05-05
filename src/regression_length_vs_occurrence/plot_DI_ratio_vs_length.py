import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression


xlsx_path = os.path.join("..", "..", "data", "schwartz2016", "SchwartzLowen_Fig3a_MdckP3.xlsx")
df = pd.read_excel(xlsx_path, skiprows=[9,10,11,12])

x = np.array(df["Length"]).reshape((-1, 1))
y = np.array(df["average"])
err = np.array(df["standard_deviation"])
segments = np.array(df["Reference"])
segments[5] = "NA"

fig, ax = plt.subplots(figsize = (10, 10))
ax.scatter(x, y, s=60, alpha=0.7, edgecolors="k")
ax.errorbar(x, y, yerr=err, fmt="o", capsize=5.0)

for i, seg in enumerate(segments):
    ax.annotate(seg, (x[i], y[i]))

model = LinearRegression().fit(x, y)
zero_model = LinearRegression(fit_intercept=False).fit(x, y)

model_x = np.linspace(0, 2341, num=100).reshape((-1, 1))

ax.plot(x, model.predict(x), label=f"linear model 1 (R²: {model.score(x, y):.2f})")
ax.plot(model_x, zero_model.predict(model_x), label=f"linear model 2 (R²: {zero_model.score(x, y):.2f})")

exp_x = x.reshape((1, -1))[0]
exp_model_x = model_x.reshape((1, -1))[0]
exp_model = np.polyfit(exp_x, np.log(y), 1)
ax.plot(exp_model_x, np.exp(exp_model[1]) * np.exp(exp_model[0] * exp_model_x), label="exponential model")

ax.set_xlabel("sequence length")
ax.set_ylabel("ratio terminal/internal")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
plt.legend(loc="upper left")
plt.suptitle("ratio of terminal to internal primers against sequence length")

fig_path = os.path.join("results", "regression_internal_vs_external.pdf")
plt.savefig(fig_path)
