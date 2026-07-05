import tkinter as tk

import sys
import os

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if root_path not in sys.path:
    sys.path.insert(0, root_path)

from App.front import GUI
from App.front.GUI import BioWorkbench
from App.front.themes.retro import config as retro_theme
from App.front.themes.modern import config as modern_theme

from App.front.themes.retro import config as retro_conf
from App.front.themes.retro import button as retro_btn

from App.front.themes.modern import config as modern_conf
from App.front.themes.modern import button as modern_btn

class AppLauncher:
    def __init__(self, root):
        self.root = root
        self.root.title("DigiLab Core Bootloader")
        self.root.geometry("400x250")

        tk.Label(self.root, text= "SELECT A SYSTEM INTERFACE", font=("Helvetica", 14)).pack(pady=20)

        tk.Button(self.root, text="LAUNCH RETRO WORKBENCH", width=25 ,command=lambda: self.boot_app(retro_conf, retro_btn)).pack(pady=10)
        tk.Button(self.root, text="LAUNCH MODERN WORKBENCH", width=25, command=lambda: self.boot_app(modern_conf, modern_btn)).pack(pady=10)

    def boot_app(self, selected_config, selected_btn):
        self.root.destroy()

        main_window = tk.Tk()
        app = BioWorkbench(main_window, selected_config, selected_btn)
        main_window.mainloop()

if __name__ == "__main__":
    launcher_root = tk.Tk()
    launcher = AppLauncher(launcher_root)
    launcher_root.mainloop()