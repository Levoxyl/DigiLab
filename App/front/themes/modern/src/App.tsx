import { useState, useEffect } from 'react'
import reactLogo from './assets/react.svg'
import viteLogo from './assets/vite.svg'
import heroImg from './assets/hero.png'
import './App.css'

import { QWebChannel } from './qwebchannel'

// Define an explicit structural shape for our Python Qt bridge
interface QtBackendBridge {
  process_dna_seq: (fastaPath: string, callback: (response: string) => void) => void;
}

function App() {
  const [backendStatus, setBackendStatus] = useState("System Offline: Initializing...")
  const [pyBackend, setPyBackend] = useState<QtBackendBridge | null>(null)

  useEffect(() => {
    // Check global scope objects exposed by the PyQt container window
    const extendedWindow = window as unknown as {
      qtWidget?: { webChannelTransport: unknown };
      QtWebChannel?: unknown;
    }

    if (extendedWindow.qtWidget) {
      new QWebChannel(extendedWindow.qtWidget.webChannelTransport, (channel) => {
        const backend = channel.objects.pyBack as QtBackendBridge
        setPyBackend(backend)
        
        // Use a functional state update to comply with react-hooks/set-state-in-effect
        setBackendStatus(() => "System Bridge Operational")
      })
    } else {
      setBackendStatus(() => "Running in standalone browser mode (Backend detached)")
    }
  }, [])

  const triggerBackendAnalysis = () => {
    if (pyBackend) {
      setBackendStatus("Sending layout request down the pipe...")
      
      pyBackend.process_dna_seq("Sample_Sequence.fasta", (response: string) => {
        setBackendStatus(`Python Result: ${response}`)
      })
    } else {
      alert("Backend bridge connection is currently unavailable.")
    }
  }

  return (
    <>
      <section id="center">
        <div className="hero">
          <img src={heroImg} className="base" width="170" height="179" alt="" />
          <img src={reactLogo} className="framework" alt="React logo" />
          <img src={viteLogo} className="vite" alt="Vite logo" />
        </div>
        <div>
          <h1>Modern Workbench Engine</h1>
          <p style={{ color: '#4ade80', fontFamily: 'monospace' }}>{backendStatus}</p>
        </div>

        <button type="button" className="counter" onClick={triggerBackendAnalysis}>
          TRIGGER PYQT BACKEND PROCESSING
        </button>
      </section>
    </>
  )
}

export default App