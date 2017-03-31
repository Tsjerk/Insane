import sys
import nose
from easy_copy import ImprovedDisplay

if __name__ == '__main__':
    nose.main(addplugins=[ImprovedDisplay()], argv=sys.argv)
