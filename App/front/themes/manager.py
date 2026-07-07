import os
import sys

from App.front.components.config import AppConfigComponent
from App.front.components.button import ButtonComponent

from App.front.themes.retro import retro_theme as retro_raw
from App.front.themes.modern import modern_theme as modern_raw

from App.front.themes.retro.layout import ThemeLayout as RetroLayout
from App.front.themes.modern.layout import ThemeLayout as ModernLayout

class ThemePack:
    """Combines configured UI components into a single, cohesive layout context."""
    def __init__(self, raw_config, raw_button, layout_class, render_engine="tkinter"):
        self.config = AppConfigComponent(raw_config)
        self.button = ButtonComponent(raw_button)
        self.layout_class = layout_class
        self.render_engine = render_engine

class ThemeManager:
    def __init__(self):
        self._registry = {
            "retro" : ThemePack(retro_raw, retro_raw, RetroLayout, render_engine="tkinter"),
            "modern" : ThemePack(modern_raw, modern_raw, ModernLayout, render_engine="chromium")
        }
        
    def get_theme(self, theme_id):
        return self._registry.get(theme_id, self._registry["retro"])
    
    def launch_theme_engine(self, theme_id, controller):
        """
        Evaluates the rendering target profile and establishes the 
        target system environment execution loop 
        """
        theme_pack = self.get_theme(theme_id)
       
        if theme_pack.render_engine == "chromium":
            from PyQt6.QtWidgets import QApplication, QMainWindow

            os.environ["QT_API"] = "pyqt6"
            os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
            os.environ["QTWEBENGINE_REMOTE_DEBUGGING"] = "9222"
            
            cleaned_args = [
                "--no-sandbox",
                "--disable-web-security",
                "--disable-gpu"
            ]

            # Combined system args with custom flags to prevent internal chromium launch failures
            app = QApplication(sys.argv + cleaned_args)
            main_window = QMainWindow()

            main_window.setWindowTitle(theme_pack.config.window_title)

            ui = theme_pack.layout_class(main_window, theme_pack, controller)
            ui.build_interface()
            main_window.showMaximized()
            sys.exit(app.exec())
        else:
            import tkinter as tk
            from App.front.GUI import BioWorkbench

            main_window = tk.Tk()
            app = BioWorkbench(main_window, theme_pack)
            main_window.mainloop()