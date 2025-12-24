#!/usr/bin/env python3
"""Normalize empirical names used in empirical overlays."""


def normalize_name_for_empirical_matching(name, entity_type='bacteria', data_source=None):
    """Normalize bacteria/drug identifiers prior to empirical lookups."""
    if name is None:
        return None

    if entity_type == 'drug':
        regions = ['north_america', 'south_america', 'europe', 'asia', 'africa', 'oceania']
        normalized = name
        for region in regions:
            prefix = f"{region}_"
            if normalized.startswith(prefix):
                normalized = normalized[len(prefix):]
                break
        return normalized

    # Bacteria names are expected to already match the canonical slugs
    return str(name).strip()