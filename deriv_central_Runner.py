import os

for i in range(1, 11):
    for prob in [25, 50, 75]:
        os.system("python deriv_centralESP.1.py " + str(prob) + " " + str(i))

        os.system("python deriv_centralFL.1.py " + str(prob) + " " + str(i))
