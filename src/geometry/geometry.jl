#! /usr/bin/julia

using GeometricTypes

# Define new abstract types for declaring 
abstract absSolidComp
abstract absVirtualComp
abstract absBoundaryComp <: absVirtualComp

