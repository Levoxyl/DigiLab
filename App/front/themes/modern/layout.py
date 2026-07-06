import os
from PyQt6.QtCore import QObject, pyqtSlot, QUrl
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWebChannel import QWebChannel

from App.back.Processing import ProcessManager

class ModernBridge(QObject):
    """ The direct messaging pipe between JavaScript and the Python backend processing """
    @pyqtSlot(str, result=str)
    def process_dna_seq(self, fasta_path):
        print(f"\n [BRIDGE SUCCESS] Python CPU recieved data from Frontend: {fasta_path}")
        result = ProcessManager.parse_dna(fasta_path)
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

        if os.path.exists(dist_path):
            self.browser.setUrl(QUrl.fromLocalFile(dist_path))
        else:
            print(f"[Vite Dist Missing] Loading modern core testing layout harness...")
            test_html = """
            <!DOCTYPE html>
            <html>
            <head>
                <style>
                    body {
                        background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
                        color: white;
                        font-family: 'Segoe UI', sans-serif;
                        display: flex;
                        flex-direction: column;
                        align-items: center;
                        justify-content: center;
                        height: 100vh;
                        margin: 0;
                        overflow: hidden;
                    }
                    h1 { font-weight: 300; letter-spacing: 2px; margin-bottom: 5px; }
                    p { color: #cbd5e1; margin-bottom: 30px; }
                    button {
                        background: #3b82f6;
                        color: white;
                        border: none;
                        padding: 12px 24px;
                        font-size: 16px;
                        font-weight: bold;
                        border-radius: 8px;
                        cursor: pointer;
                        box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1);
                        transition: all 0.2s;
                    }
                    button:hover { background: #2563eb; transform: translateY(-1px); }
                    #output { margin-top: 20px; color: #4ade80; font-family: monospace; font-size: 14px; }
                </style>
                <script src="qrc:///qtwebchannel/qwebchannel.js"></script>
            </head>
            <body>
                <h1>CHROMIUM CORE ACTIVE</h1>
                <p>Modern Engine Workspace Operational (React/TS/Three.js Pipeline Setup)</p>
                <button id="testBtn">TRIGGER PYTHON CPU TEST MESSAGE</button>
                <div id="output"></div>

                <script>
                    document.getElementById('testBtn').addEventListener('click', () => {
                        if (window.qtWidget || window.QtWebChannel) {
                            new QtWebChannel(qtWidget.webChannelTransport, function(channel) {
                                let backend = channel.objects.pyBack;
                                document.getElementById('output').innerText = "Sending payload to Python...";
                                
                                // Call our backend method securely
                                backend.process_dna_seq("Test_Sample_Sequence.fasta", function(response) {
                                    document.getElementById('output').innerText = "Python Response: " + response;
                                });
                            });
                        } else {
                            document.getElementById('output').innerText = "Error: QWebChannel channel transport missing.";
                        }
                    });
                </script>
            </body>
            </html>
            """
            self.browser.setHtml(test_html)

        self.root.setCentralWidget(self.browser)
        self.browser.show()


        # self.root.title(self.theme.window_title)
        # self.root.geometry(self.theme.window_size)

        # browser_frame = tk.Frame(self.root)
        # browser_frame.pack(fill="both", expand=True)