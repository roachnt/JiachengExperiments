import os

for prob in [25, 50, 75]:
    os.system("python J0LongESP.py " + str(prob))
    os.system("python J0LongESP.1.py " + str(prob))
    os.system("python J0LongESP.2.py " + str(prob))
    os.system("python J0LongESP.3.py " + str(prob))
    os.system("python J0LongESP.4.py " + str(prob))

    os.system("python J0LongFL.py " + str(prob))
    os.system("python J0LongFL.1.py " + str(prob))
    os.system("python J0LongFL.2.py " + str(prob))
    os.system("python J0LongFL.3.py " + str(prob))
    os.system("python J0LongFL.4.py " + str(prob))
