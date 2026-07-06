class AppConfigComponent:
    def __init__(self, raw_config_module):
        self.window_title = getattr(raw_config_module, "WINDOW_TITLE", "Genomic Analyzer")
        self.window_size = getattr(raw_config_module, "WINDOW_SIZE", "800x600")
        self.bg_color = getattr(raw_config_module, "BG_COLOR", "#ffffff")

        self.header_style = getattr(raw_config_module, "HEADER_STYLE", {})
        self.log_style = getattr(raw_config_module, "LOG_STYLE", {})