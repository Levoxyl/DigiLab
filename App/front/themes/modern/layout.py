import os
from PyQt6.QtCore import QObject, pyqtSlot, QUrl
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWebChannel import QWebChannel

from App.back.Processing import Processing

class ModernBridge(QObject):
    """ The direct messaging pipe between JavaScript and the Python backend processing """
    @pyqtSlot(str, result=str)
    def process_dna_seq(self, fasta_path):
        result = Processing.parse_dna(fasta_path)
        return result

class ThemeLayout:
    """Implements an embedded Chromium browser view frame."""
    def __init__(self, root, theme_pack, controller):
        self.root = root
        self.theme = theme_pack.config
        self.controller = controller

    def build_interface(self):
        self.browser = QWebEngineView(self.root)
        self.channel = QWebChannel()
        self.bridge = ModernBridge()

        self.channel.registerObject("pyBack", self.bridge)
        self.browser.page().setWebChannel(self.channel)

        dist_path = os.path.abspath(os.path.join(os.path.dirname(__file__),"dist", "index.html"))
        self.browser.setUrl(QUrl.fromLocalFile(dist_path))

        self.browser.show()


        self.root.title(self.theme.window_title)
        self.root.geometry(self.theme.window_size)

        browser_frame = tk.Frame(self.root)
        browser_frame.pack(fill="both", expand=True)