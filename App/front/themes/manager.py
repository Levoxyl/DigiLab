from App.front.components.config import AppConfigComponent
from App.front.components.button import ButtonComponent

from App.front.themes.retro import config as retro_raw_conf
from App.front.themes.retro import button as retro_raw_btn
from App.front.themes.retro.layout import ThemeLayout as Retrolayout 

from App.front.themes.modern import config as modern_raw_conf
from App.front.themes.modern import button as modern_raw_btn
from App.front.themes.modern.layout import ThemeLayout as ModernLayout

class ThemePack:
    """ Combines configurated UI components into a single ,cohesive layout context. """
    def __init__(self,raw_config, raw_button, layout_class, render_engine="tkinter"):
        self.config = AppConfigComponent(raw_config)
        self.button = ButtonComponent(raw_button)
        self.render_engine = render_engine
        self.layout_class = layout_class

class ThemeManager:
    def __init__(self):
        self._registry = {
            "retro" : ThemePack(retro_raw_conf, retro_raw_btn, Retrolayout, render_engine="tkinter"),
            "modern" : ThemePack(modern_raw_conf, modern_raw_btn, ModernLayout, render_engine="chromium")
        }
    def get_theme(self, theme_id):
        return self._registry.get(theme_id, self._registry["retro"])