import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import type { DecisionResult } from "@/lib/api"

interface DecisionSummaryProps {
  result: DecisionResult
}

const DECISION_STYLES: Record<string, { bg: string; text: string }> = {
  TREAT: { bg: "bg-green-50 border-green-200", text: "text-green-700" },
  "DO NOT TREAT": { bg: "bg-slate-50 border-slate-200", text: "text-slate-700" },
  "REQUEST EXAMS": { bg: "bg-amber-50 border-amber-200", text: "text-amber-700" },
}

export default function DecisionSummary({ result }: DecisionSummaryProps) {
  const style = DECISION_STYLES[result.decision] ?? {
    bg: "bg-secondary border-border",
    text: "text-foreground",
  }
  const { mu, nu, pi } = result.consensus

  return (
    <Card className="border-border shadow-sm hover:shadow-md transition-shadow">
      <CardHeader>
        <CardTitle>Clinical Recommendation</CardTitle>
        <CardDescription>Analysis based on AI assessments</CardDescription>
      </CardHeader>

      <CardContent className="space-y-6">
        <div className={`p-6 rounded-lg border-2 ${style.bg}`}>
          <p className={`text-xs uppercase font-bold tracking-widest mb-2 ${style.text}`}>
            Final Recommendation
          </p>
          <p className={`text-3xl font-bold ${style.text}`}>{result.decision}</p>
        </div>

        <div className="pt-4 border-t border-border">
          <p className="text-sm font-semibold text-foreground mb-3">Clinical Reasoning</p>
          <div className="bg-primary/5 border border-primary/10 rounded-lg p-4">
            <p className="text-sm leading-relaxed text-foreground">
              μ⊙ = {mu.toFixed(2)}, ν⊙ = {nu.toFixed(2)}, π⊙ = {pi.toFixed(2)}. (μ⊙ =
              agreement, ν⊙ = conflict, π⊙ = hesitation). Final decision:{" "}
              {result.decision}.
            </p>
          </div>
        </div>
      </CardContent>
    </Card>
  )
}
