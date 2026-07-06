import os

from App.front.components.config import AppConfigComponent
from App.front.components.button import ButtonComponent

# Pull raw theme parameters from their new home inside the themes/ directory
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
        target system enviroment execution loop 
        """
        theme_pack = self.get_theme(theme_id)
       
        if theme_pack.render_engine == "chromium":
                import sys
                import os
                from PyQt6.QtWidgets import QApplication, QMainWindow

               # Force pure software fallback rendering without GPU hooks
                os.environ["QTWEBENGINE_DISABLE_SANDBOX"] = "1"
                os.environ["QT_X11_NO_MITSHM"] = "1"
                os.environ["QT_DISABLE_WEBGL_HARDWARE_ACCELERATION"] = "1"
                os.environ["QT_QUICK_BACKEND"] = "software"
                os.environ["QSG_RHI_BACKEND"] = "softwarecontext" 
                os.environ["QT_API"] = "pyqt6"

                # Keep your existing fallback args to catch background loops
                sys.argv.append("--no-sandbox")
                sys.argv.append("--disable-gpu")
                sys.argv.append("--disable-software-rasterizer")
                sys.argv.append("--ignore-gpu-blocklist")

                app = QApplication(sys.argv)
                main_window = QMainWindow()

                main_window.setWindowTitle(theme_pack.config.window_title)

                dims = theme_pack.config.window_size.split("x")
                if len(dims) == 2:
                    main_window.resize(int(dims[0]), int(dims[1]))

                ui = theme_pack.layout_class(main_window, theme_pack, controller)
                ui.build_interface()

                main_window.show()
                sys.exit(app.exec())
        else:
            import tkinter as tk
            from App.front.GUI import BioWorkbench

            main_window = tk.Tk()
            app = BioWorkbench(main_window, theme_pack)
            main_window.mainloop()
