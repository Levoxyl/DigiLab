import React, { useState, useEffect } from 'react';

interface AbstractRiverBackgroundProps {
  children?: React.ReactNode;
}

const AbstractRiverBackground: React.FC<AbstractRiverBackgroundProps> = ({ children }) => {
  const [renderStage, setRenderStage] = useState<number>(0);

  useEffect(() => {
    const tick1 = setTimeout(() => setRenderStage(1), 30);
    const tick2 = setTimeout(() => setRenderStage(2), 150);
    const tick3 = setTimeout(() => setRenderStage(3), 400);

    return () => {
      clearTimeout(tick1);
      clearTimeout(tick2);
      clearTimeout(tick3);
    };
  }, []);

  return (
    <div className="freeform-gradient-wrapper" data-stage={renderStage}>

      {/* CENTRALIZED HARDWARE LAYER CACHE */}
      <div className="gpu-canvas-cache-container">
        <svg
          className="abstract-fluid-matrix"
          xmlns="http://www.w3.org/2000/svg"
          viewBox="0 0 1920 1080"
          preserveAspectRatio="xMidYMid slice"
        >
          <rect width="100%" height="100%" fill="#b2cfc6" />

          {/* MAIN COLOR STREAM */}
          <path
            className="vector-blob blob-navy stream-primary"
            fill="#223861"
            style={{ mixBlendMode: 'multiply' }}
            opacity="0.80"
            d="M960,540 C1250,200 1600,300 1750,550 C1900,800 1450,1050 960,950 C470,850 300,750 250,550 C200,350 670,880 960,540 Z"
          />
          <path
            className="vector-blob blob-teal-dark stream-primary"
            fill="#1a5c59"
            style={{ mixBlendMode: 'multiply' }}
            opacity="0.70"
            d="M800,400 C1150,150 1450,300 1550,550 C1650,800 1250,950 850,850 C450,750 350,650 300,450 C250,250 450,650 800,400 Z"
          />
          <path
            className="vector-blob blob-teal-light stream-primary"
            fill="#2a7a77"
            style={{ mixBlendMode: 'screen' }}
            opacity="0.60"
            d="M1100,600 C1350,350 1650,400 1700,650 C1750,900 1400,1000 1100,900 C800,800 650,700 650,550 C650,400 850,850 1100,600 Z"
          />

          {/* LOWER LATERAL BALANCING SYSTEM */}
          <g transform="translate(-400, 300) scale(0.9)">
            <path
              className="vector-blob blob-teal-light stream-tertiary"
              fill="#2a7a77"
              style={{ mixBlendMode: 'screen' }}
              opacity="0.50"
              d="M1100,600 C1350,350 1650,400 1700,650 C1750,900 1400,1000 1100,900 C800,800 650,700 650,550 C650,400 850,850 1100,600 Z"
            />
          </g>
        </svg>

        {/* RE-INTEGRATED HIGH-DPI FILM GRAIN SYSTEM */}
        <div className="cinematic-grain-overlay" />
      </div>

      {/* ISOLATED INTERFACE WRAPPER */}
      <div className="background-content-container">
        {children}
      </div>
    </div>
  );
};

export default AbstractRiverBackground;