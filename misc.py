import sys

def msgError(str, numError = 1):
    print('*--------------------------------------------------------')
    print(str)
    print('*--------------------------------------------------------')
    sys.exit(numError)


def exFloatStrToInt(strFloat):
    try:
        return float(strFloat)
    except:
        # If it doesn't work check if it is a float that uses D instead of E...
        return float(strFloat.lower().replace('d', 'e'))
