"use client"

import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { hesitation, isValid, type Opinion } from "@/lib/api"

interface AiOpinionsProps {
  opinions: Opinion[]
  onChange: (opinions: Opinion[]) => void
}

/**
 * Defined at module scope on purpose: declaring it inside the parent would
 * create a new component type on every render, remounting the input and
 * dropping keyboard focus while dragging.
 */
function SliderInput({
  label,
  symbol,
  value,
  onValueChange,
}: {
  label: string
  symbol: string
  value: number
  onValueChange: (value: number) => void
}) {
  return (
    <div className="space-y-2.5">
      <div className="flex items-center justify-between">
        <label className="text-sm font-medium text-foreground">
          {label} <span className="text-primary font-semibold">{symbol}</span>
        </label>
        <span className="text-sm font-semibold text-primary bg-primary/10 px-2.5 py-1 rounded">
          {value.toFixed(2)}
        </span>
      </div>
      <input
        type="range"
        min="0"
        max="1"
        step="0.01"
        value={value}
        onChange={(e) => onValueChange(Number.parseFloat(e.target.value))}
        className="w-full h-2.5 bg-secondary rounded-full appearance-none cursor-pointer accent-primary"
      />
    </div>
  )
}

export default function AiOpinions({ opinions, onChange }: AiOpinionsProps) {
  const update = (id: string, patch: Partial<Opinion>) =>
    onChange(opinions.map((o) => (o.id === id ? { ...o, ...patch } : o)))

  return (
    <Card className="border-border shadow-sm hover:shadow-md transition-shadow">
      <CardHeader>
        <CardTitle>AI Opinions</CardTitle>
        <CardDescription>Set belief measures from each AI</CardDescription>
      </CardHeader>

      <CardContent className="space-y-8">
        {opinions.map((opinion, index) => {
          const invalid = !isValid(opinion)
          const pi = hesitation(opinion)

          return (
            <div
              key={opinion.id}
              className={
                index < opinions.length - 1
                  ? "space-y-4 pb-6 border-b border-border/50"
                  : "space-y-4"
              }
            >
              <div className="flex items-center gap-2 mb-4">
                <div className="w-8 h-8 shrink-0 rounded-full bg-primary/10 border border-primary/20 flex items-center justify-center">
                  <span className="text-sm font-bold text-primary">{index + 1}</span>
                </div>
                <h3 className="font-semibold text-foreground">{opinion.label}</h3>
              </div>

              <SliderInput
                label="Belief in Disease"
                symbol={`μ${index + 1}`}
                value={opinion.mu}
                onValueChange={(mu) => update(opinion.id, { mu })}
              />
              <SliderInput
                label="Belief in NO Disease"
                symbol={`ν${index + 1}`}
                value={opinion.nu}
                onValueChange={(nu) => update(opinion.id, { nu })}
              />

              {/* Hesitation band: the mass left over up to 1 */}
              <div
                className={`p-3 rounded-lg border-2 transition-colors ${
                  invalid
                    ? "bg-destructive/10 border-destructive/30 text-destructive"
                    : "bg-primary/5 border-primary/20 text-foreground"
                }`}
              >
                <p className="text-xs font-medium">
                  {invalid
                    ? "Invalid: μ + ν cannot exceed 1"
                    : `Hesitation π${index + 1} = ${pi.toFixed(2)}`}
                </p>
              </div>
            </div>
          )
        })}
      </CardContent>
    </Card>
  )
}
