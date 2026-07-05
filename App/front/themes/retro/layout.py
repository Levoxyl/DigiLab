import tkinter as tk

class ThemeLayout:
    """ Implements a pure Tkinter desktop grid rendering for retro theme """
    def __init__(self, root, theme_pack, controller):
        self.root = root
        self.theme = theme_pack.config
        self.button_theme = theme_pack.button
        self.controller = controller

    def build_interface(self):
        
        # Base Framework Generation
        self.root.title(self.theme.window_title)
        self.root.geometry(self.theme.window_size)
        self.root.configure(bg=self.theme.bg_color)

        # Header Frame
        header = tk.Label(self.root, text="LAB_TERMINAL: DNA_DIGEST_SYSTEM", **self.theme.header_style)
        header.pack(fill="x", pady=5)

        # 3. Input Panels Section
        input_frame = tk.Frame(self.root, bg=self.theme.bg_color, bd=2, relief="groove")
        input_frame.pack(fill="x", padx=10, pady=10)

        # Enforce a custom width just for these two structural buttons so they stay identical
        tk.Button(input_frame, text="FILE...", command=self.controller.select_file, **{**self.button_theme.style, "width": 10}).grid(row=0, column=0, padx=8, pady=5)
        tk.Entry(input_frame, textvariable=self.controller.fasta_path, width=60, bg="white").grid(row=0, column=1, padx=5)

        tk.Button(input_frame, text="DIR...", command=self.controller.select_folder, **{**self.button_theme.style, "width": 10}).grid(row=1, column=0, padx=8, pady=5)
        tk.Entry(input_frame, textvariable=self.controller.output_dir, width=60, bg="white").grid(row=1, column=1, padx=5)

        # Log Output Feed Console
        self.controller.log_box = tk.Text(self.root, height=15, **self.theme.log_style)
        self.controller.log_box.pack(padx=10, pady=5, fill="both")

        # Bottom Operation Controls Button Row
        btn_f = tk.Frame(self.root, bg=self.theme.bg_color)
        btn_f.pack(pady=10)

        # Use ipadx (internal horizontal padding) and ipady (internal vertical padding) to make the buttons large, clean, and unclipped
        tk.Button(btn_f, text="[ RUN_TRANSLATE ]", command=self.controller.do_translate, **self.button_theme.style).pack(side="left", padx=4, ipadx=12, ipady=4)
        tk.Button(btn_f, text="[ RUN_DIGEST ]", command=self.controller.do_digest, **self.button_theme.style).pack(side="left", padx=4, ipadx=12, ipady=4)
        tk.Button(btn_f, text="[ CLEAR ]", command=self.controller.clear_log, **self.button_theme.style).pack(side="left", padx=4, ipadx=12, ipady=4)
        tk.Button(btn_f, text="[ ANALYZE_CHEMISTRY ]", command=self.controller.analyze_chemistry, **self.button_theme.style).pack(side="left", padx=4, ipadx=12, ipady=4)