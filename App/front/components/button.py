class ButtonComponent:
    def __init__(self, raw_button_module):
        self.style = getattr(raw_button_module, "STYLE", {})