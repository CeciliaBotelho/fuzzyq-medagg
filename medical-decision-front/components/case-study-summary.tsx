import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { accuracy, hesitation, score, type Opinion } from "@/lib/api"

interface CaseStudySummaryProps {
  opinions: Opinion[]
}

export default function CaseStudySummary({ opinions }: CaseStudySummaryProps) {
  const ranked = [...opinions].sort((a, b) => score(a) - score(b))
  const tied =
    ranked.length === 2 && Math.abs(score(ranked[0]) - score(ranked[1])) < 1e-12

  return (
    <Card className="border-border shadow-sm">
      <CardHeader>
        <CardTitle className="text-lg">Case Study Summary</CardTitle>
        <CardDescription>Fixed inputs for the paper example</CardDescription>
      </CardHeader>

      <CardContent className="space-y-4">
        {/* The intuitionistic fuzzy values themselves */}
        <div className="bg-secondary/40 border border-border/50 rounded-lg p-4">
          <p className="text-sm font-semibold text-foreground mb-3">
            AI opinions (Atanassov IFVs)
          </p>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-border text-muted-foreground">
                  <th className="text-left py-2 pr-2 font-medium">Opinion</th>
                  <th className="text-right py-2 px-2 font-medium">μ</th>
                  <th className="text-right py-2 px-2 font-medium">ν</th>
                  <th className="text-right py-2 pl-2 font-medium">π = 1 − μ − ν</th>
                </tr>
              </thead>
              <tbody>
                {opinions.map((o) => (
                  <tr key={o.id} className="border-b border-border/40 last:border-0">
                    <td className="py-2 pr-2 text-foreground whitespace-nowrap">{o.label}</td>
                    <td className="py-2 px-2 text-right tabular-nums">{o.mu.toFixed(2)}</td>
                    <td className="py-2 px-2 text-right tabular-nums">{o.nu.toFixed(2)}</td>
                    <td className="py-2 pl-2 text-right tabular-nums">
                      {hesitation(o).toFixed(2)}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>

        {/* Xu-Yager total order, kept deliberately secondary to the IFVs above */}
        <div className="pt-1">
          <p className="text-xs font-medium text-muted-foreground mb-1.5">
            Admissible-order comparison
          </p>
          <p className="text-xs text-muted-foreground leading-relaxed">
            Xu&ndash;Yager total order: compare s = &mu; &minus; &nu; first, then
            h = &mu; + &nu; on a tie.{" "}
            {opinions.map((o, k) => (
              <span key={o.id} className="tabular-nums">
                {k > 0 && " · "}
                {o.label}: s = {score(o).toFixed(2)}, h = {accuracy(o).toFixed(2)}
              </span>
            ))}
            {ranked.length === 2 && (
              <>
                {" "}&rarr;{" "}
                <span className="font-medium text-foreground">
                  {ranked[0].label} ≺ {ranked[1].label}
                </span>
                {tied ? ", decided by h." : ", the scores do not tie."}
              </>
            )}
          </p>
        </div>
      </CardContent>
    </Card>
  )
}
