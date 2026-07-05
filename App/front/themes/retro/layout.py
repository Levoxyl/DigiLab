import tkinter as tk

class ThemeLayout:
    """ Implements a pure Tkinter desktop grid rendering for retro theme """
    def __init__ (self,root,theme_pack,controller):
        self.root = root
        self.theme = theme_pack.config
        self.button_theme = theme_pack.button
        self.controller = controller

    def build_interface(self):
        self.root.title(self.theme.window_title)
        self.root.geometry(self.theme.window_size)
        self.root.configure(bg=self.theme.bg_color)

        header = tk.Label(self.root, text="LAB_TERMINAL: DNA_DIGEST_SYSTEM", **self.theme.header_style)
        header.pack(fill="x", pady=5)