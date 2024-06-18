import sys
from GUI import *
import time


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = HeliostatUI()
    window.show()
    sys.exit(app.exec_())
