import save
import numpy as np
import matplotlib.pyplot as plt

mu = 0.1e31

theta_obs = 60
beta_mu = 20
mc2 = 60
a_portion = 0.8
phi_0 = 0

theta_obs = 80
beta_mu = 20
mc2 = 60
a_portion = 0.2
phi_0 = 40

L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

L_total = L_total / np.max(L_total)


def get_delta(left, right):
    return np.sqrt(np.sum((left - right) ** 2)) # + np.max(np.abs(left - right))


deltas = np.zeros(L_total.shape[0])
deltas[0] = float('inf')
for i in range(1, L_total.shape[0]):
    # отражаем куски относительно фазы разреза. и дополняем до полного профиля
    # правый едет направо, левый - налево от места разреза.
    # если доходим до конца то прыгыаем на другой конец и идем до места разреза
    left = np.concatenate((L_total[:i + 1][::-1], L_total[i + 1:][::-1]))
    right = np.concatenate([L_total[i:], L_total[:i]])

    deltas[i] = get_delta(left, right)

ids = np.argsort(deltas)
print(deltas)
print(ids)

i = ids[0]

left = np.concatenate((L_total[:i + 1][::-1], L_total[i + 1:][::-1]))
right = np.concatenate([L_total[i:], L_total[:i]])
print(get_delta(left, right))
# left = np.concatenate([L_total[i:], L_total[:i]])
# right = np.concatenate([L_total[i:][::-1], L_total[:i]][::-1])[::-1]

plt.plot(L_total, label='original')
plt.plot(left, label='left')
plt.plot(right, label='right')
plt.legend()
plt.show()

# -------------------------


# L_total = np.linspace(1, 45, 45)
#
# deltas = np.zeros(L_total.shape[0])
# deltas[0] = float('inf')
# for i in range(1, L_total.shape[0]):
#     left = np.concatenate((L_total[:i + 1][::-1], L_total[i + 1:][::-1]))
#     right = np.concatenate([L_total[i:], L_total[:i]])
#
#     deltas[i] = get_delta(left, right)
#
# ids = np.argsort(deltas)
# print(deltas)
# print(ids)
#
# i = 11
# print(i)
# # roll = вправо
# print(np.roll(L_total, i))
#
# # print()
# left = np.concatenate((L_total[:i + 1][::-1], L_total[i + 1:][::-1]))
# right = np.concatenate([L_total[i:], L_total[:i]])
# print(get_delta(left, right))
# # left = np.concatenate([L_total[i:], L_total[:i]])
# # right = np.concatenate([L_total[i:][::-1], L_total[:i]][::-1])[::-1]
#
# plt.plot(L_total)
# plt.plot(left, label='left')
# plt.plot(right, label='right')
# plt.legend()
# plt.show()
#
# print(left)
# print(right)
# # print(np.stack(, axis=0))
