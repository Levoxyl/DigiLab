import React from 'react';

const LoginMenu: React.FC = () => {
    return (
        <div className = "login-bg-container">
            {/* Freeroam grid */}
            <div className="mesh-blob blob1"></div>
            <div className="mesh-blob blob2"></div>
            <div className="mesh-blob blob3"></div>

            {/* Centered Menu  Panel */}
            <div className="menu-card">
                <header className="menu-header">
                    <h1> GENOMIC ANALYZER </h1>
                    <p className="menu-subtitle">SYSTEM MANEGMENT INTERFACE</p>
                </header>

                {/* Sharp cornered gradient buttons */}
                <div className="button-group">
                    <button type="button"  className="sharp-gradient-btn">
                        LAUNCH RETRO WORKBENCH
                    </button>
                    <button type="button"  className="sharp-gradient-btn">
                        RUN SEQUENCE ANALYZER
                    </button>
                    <button type="button"  className="sharp-gradient-btn">
                        VIEW CONNECTION LOGS
                    </button>
                </div>
            </div>
        </div>
    );
};

export default LoginMenu;