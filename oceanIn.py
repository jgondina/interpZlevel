import configparser
import textwrap
from misc import *


class OceanIn:
    def __init__(self, fileName):
        strOceanIn = ''
        # Read an ocean.in file removing the leading spaces.
        fileOceanin = open(fileName)
        strOceanInOriginal = textwrap.dedent(fileOceanin.read())
        for line in strOceanInOriginal.splitlines():
            strOceanIn += line.lstrip() + '\n'
        # Fakes a section
        strOceanIn = '[sec]\n' + strOceanIn

        parser = configparser.ConfigParser(comment_prefixes=('!'), inline_comment_prefixes=('!'), delimiters=('==', '='))

        parser.read_string(strOceanIn)

        self.theta_s = exFloatStrToInt(parser['sec']['theta_s'])
        self.theta_b = exFloatStrToInt(parser['sec']['theta_b'])
        self.tcline  = exFloatStrToInt(parser['sec']['tcline'])

        self.N           = int(parser['sec']['N'])

        self.Vtransform  = int(parser['sec']['Vtransform'])
        self.Vstretching = int(parser['sec']['Vstretching'])

        pass
