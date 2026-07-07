import React from 'react';
import AbstractRiverBackground from '../styles/uiComponents/AbstractRiverBackground';

interface LoginMenuProps {
  backendStatus: string;
  onTriggerAnalysis: () => void;
}

const LoginMenu: React.FC<LoginMenuProps> = ({ backendStatus, onTriggerAnalysis }) => {
  return (
    <AbstractRiverBackground>
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
    </AbstractRiverBackground>
  );
};

export default LoginMenu;