import sys
import nose
from color import ColorVerbose

if __name__ == '__main__':
    nose.main(addplugins=[ColorVerbose()], argv=sys.argv)
