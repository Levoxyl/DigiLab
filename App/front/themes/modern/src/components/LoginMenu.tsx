import React from 'react';

interface LoginMenuProps {
  backendStatus: string;
  onTriggerAnalysis: () => void;
}

const LoginMenu: React.FC<LoginMenuProps> = ({ backendStatus, onTriggerAnalysis }) => {
  return (
    <div className="login-bg-container">
      {/* Film Grain Layer Overlays Everything */}
      <div className="film-grain-overlay"></div>

      {/* Freeform Smooth Mesh Gradient Spheres */}
      <div className="mesh-blob blob1"></div>
      <div className="mesh-blob blob2"></div>
      <div className="mesh-blob blob3"></div>

      {/* Clean Left-Aligned Core Layout - No Backdrop Wrapper box */}
      <div className="dashboard-layout-engine">
        <header className="menu-header">
          <h1>GENOMIC ANALYZER</h1>
          <p className="menu-subtitle">{backendStatus}</p>
        </header>

        <div className="button-group">
          <button type="button" className="sharp-gradient-btn">
            LAUNCH RETRO WORKBENCH
          </button>
          
          <button type="button" className="sharp-gradient-btn" onClick={onTriggerAnalysis}>
            RUN SEQUENCE ANALYZER
          </button>
          
          <button type="button" className="sharp-gradient-btn">
            VIEW CONNECTION LOGS
          </button>
        </div>
      </div>
    </div>
  );
};

export default LoginMenu;