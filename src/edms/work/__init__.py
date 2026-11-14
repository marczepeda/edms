'''
Subpackage: work
Author: Marc Zepeda
Created: 2026-11-12
Description: Aided experiment planning and process tracking.

Turns a protocol template plus a short run spec into a concrete, calculated
protocol, then keeps a ledger of what actually happened as you work through it.

Usage:
- template: load and validate protocol templates
- spec: parse run spec files
- plan: template + spec -> resolved protocol (calls into edms.bio)
- run: persistent run state (status, notes, deviations)
- render: resolved protocol + run state -> Markdown / HTML
- plates: culture vessel formats and plate-aware transfection planning
'''
from . import template, spec, plan, run, render, plates

__all__ = ['template', 'spec', 'plan', 'run', 'render', 'plates']
