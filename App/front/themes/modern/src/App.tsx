import { useState, useEffect } from 'react'
import './App.css'
import { QWebChannel } from './qwebchannel'
import LoginMenu from './components/LoginMenu'

interface QtBackendBridge {
  process_dna_seq: (fastaPath: string, callback: (response: string) => void) => void;
}

function App() {
  // Fix 1: Make the standalone text the default state to avoid synchronous updates inside useEffect
  const [backendStatus, setBackendStatus] = useState("Running in standalone browser mode (Backend detached)")
  const [pyBackend, setPyBackend] = useState<QtBackendBridge | null>(null)

  useEffect(() => {
    const extendedWindow = window as unknown as {
      qtWidget?: { webChannelTransport: unknown };
      QtWebChannel?: unknown;
    }

    if (extendedWindow.qtWidget) {
      new QWebChannel(extendedWindow.qtWidget.webChannelTransport, (channel) => {
        const backend = channel.objects.pyBack as QtBackendBridge
        setPyBackend(backend)
        setBackendStatus(() => "System Bridge Operational")
      })
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

  // Fix 2: Pass down the tracking values to clear the unused variable linter errors
  return (
    <LoginMenu 
      backendStatus={backendStatus} 
      onTriggerAnalysis={triggerBackendAnalysis}
    />
  )
}

export default App