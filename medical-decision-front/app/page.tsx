"use client"

import { useState } from "react"
import Header from "@/components/header"
import AiOpinions from "@/components/ai-opinions"
import DecisionSummary from "@/components/decision-summary"
import CaseStudySummary from "@/components/case-study-summary"
import { Button } from "@/components/ui/button"
import {
  isValid,
  requestDecision,
  type CombineMode,
  type DecisionResult,
  type Opinion,
} from "@/lib/api"

// Fixed inputs of the paper example.
const INITIAL_OPINIONS: Opinion[] = [
  { id: "ai-1", label: "AI 1", mu: 0.85, nu: 0.1 },
  { id: "ai-2", label: "AI 2", mu: 0.8, nu: 0.02 },
]

export default function Page() {
  const [opinions, setOpinions] = useState<Opinion[]>(INITIAL_OPINIONS)
  const [combine, setCombine] = useState<CombineMode>("uniform")
  const [result, setResult] = useState<DecisionResult | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(false)

  const allValid = opinions.every(isValid)

  const handleCompute = async () => {
    setLoading(true)
    setError(null)
    try {
      // Closed form: the figure shows the exact degrees of Proposition 2,
      // free of the sampling noise that would make the last decimal wobble.
      setResult(await requestDecision(opinions, 20000, undefined, combine, true))
    } catch (err) {
      setResult(null)
      setError(
        err instanceof Error
          ? err.message
          : "Could not reach the backend. Make sure it is running.",
      )
    } finally {
      setLoading(false)
    }
  }

  return (
    <div className="min-h-screen bg-background">
      <Header />

      <main className="flex-1">
        <div className="container mx-auto px-4 py-10 md:py-14">
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
            {/* Left column - input controls */}
            <div className="lg:col-span-1 space-y-6">
              <AiOpinions opinions={opinions} onChange={setOpinions} />

              <Button
                onClick={handleCompute}
                disabled={!allValid || loading}
                className="w-full h-11 text-base font-semibold rounded-lg"
                size="lg"
              >
                {loading ? "Analyzing…" : "Analyze Assessment"}
              </Button>

              {!allValid && (
                <p className="text-sm text-destructive text-center">
                  Fix the highlighted opinions before analyzing.
                </p>
              )}
            </div>

            {/* Right column - results and visualizations */}
            <div className="lg:col-span-2 space-y-6">
              {error && (
                <div className="bg-destructive/10 border-2 border-destructive/30 rounded-lg p-4">
                  <p className="text-sm font-semibold text-destructive mb-1">Analysis failed</p>
                  <p className="text-sm text-destructive/90">{error}</p>
                </div>
              )}

              {result ? (
                <DecisionSummary result={result} />
              ) : (
                !error && (
                  <div className="bg-card border border-border rounded-lg p-12 text-center">
                    <div className="space-y-2">
                      <p className="text-muted-foreground text-base">
                        Set the diagnostic sources and click
                      </p>
                      <p className="font-semibold text-foreground">&quot;Analyze Assessment&quot;</p>
                      <p className="text-muted-foreground text-sm">
                        to generate a clinical recommendation
                      </p>
                    </div>
                  </div>
                )
              )}

              <CaseStudySummary opinions={opinions} />
            </div>
          </div>
        </div>
      </main>
    </div>
  )
}
