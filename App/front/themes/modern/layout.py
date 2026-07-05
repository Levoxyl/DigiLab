import tkinter as tk

class ThemeLayout:
    """Implements an embedded Chromium browser view frame."""
    def __init__(self, root, theme_pack, controller):
        self.root = root
        self.theme = theme_pack.config
        self.controller = controller

    def build_interface(self):
        self.root.title(self.theme.window_title)
        self.root.geometry(self.theme.window_size)

        browser_frame = tk.Frame(self.root)
        browser_frame.pack(fill="both", expand=True)