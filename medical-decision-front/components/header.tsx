export default function Header() {
  return (
    <header className="bg-gradient-to-br from-primary/5 to-accent/5 border-b border-border">
      <div className="container mx-auto px-4 py-12 md:py-16">
        <div className="max-w-2xl">
          <div className="inline-block mb-3 px-3 py-1 rounded-full bg-primary/10 border border-primary/20">
            <p className="text-xs font-semibold text-primary uppercase tracking-wider">AI Opinion Evaluation</p>
          </div>
          <h1 className="text-4xl md:text-5xl font-bold text-foreground mb-3 leading-tight">
            AI Opinion Assessment System
          </h1>
          <p className="text-base md:text-lg text-muted-foreground max-w-xl">
            A system for analyzing and comparing AI-generated opinions using logical evaluation and advanced interpretability methods          </p>
        </div>
      </div>
    </header>
  )
}
