import os
import sys

if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.argv.append(os.path.dirname(__file__))
    import nose
    nose.main()
