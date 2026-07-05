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