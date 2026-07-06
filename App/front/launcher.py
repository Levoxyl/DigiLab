import tkinter as tk

import sys
import os

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if root_path not in sys.path:
    sys.path.insert(0, root_path)

from App.front import GUI
from App.front.GUI import BioWorkbench

from App.front.themes.manager import ThemeManager

class AppLauncher:
    def __init__(self, root):
        self.root = root
        self.theme_manager = ThemeManager()

        tk.Label(self.root, text= "SELECT A SYSTEM INTERFACE", font=("Helvetica", 14)).pack(pady=20)

        tk.Button(self.root, text="LAUNCH RETRO WORKBENCH",
                   width=25 ,command=lambda: self.boot_app("retro")).pack(pady=10)
        tk.Button(self.root, text="LAUNCH MODERN WORKBENCH",
                   width=25, command=lambda: self.boot_app("modern")).pack(pady=10)

    def boot_app(self, theme_id):
        selected_theme_pack = self.theme_manager.get_theme(theme_id)
        self.root.destroy()

        self.theme_manager.launch_theme_engine(theme_id, controller=None)

        main_window = tk.Tk()
        app = BioWorkbench(main_window, selected_theme_pack)
        main_window.mainloop()

if __name__ == "__main__":
    launcher_root = tk.Tk()
    launcher = AppLauncher(launcher_root)
    launcher_root.mainloop()