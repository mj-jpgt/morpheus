"""Backward-compatible entry point for the complete contract-aware task suite.

The previous script selected whichever representation happened to be present
and called a seen-label mapper "zero-shot".  Keeping this wrapper ensures old
automation uses the strict comprehensive evaluator instead.
"""
from .comprehensive_evaluation import main


if __name__ == "__main__":
    main()
