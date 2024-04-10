module DebugAsserts
# Thanks are due to c42f for this!

export debug_asserts, @dassert

function debug_asserts(m::Module, dodebug::Bool)
    m.eval(:(_debug_enabled() = $dodebug))
    nothing
end

# Like @assert, but only active in debug mode
macro dassert(exs...)
    if !isdefined(__module__, :_debug_enabled)
        eval(__module__, :(_debug_enabled() = false))
    end
    quote
        if $__module__._debug_enabled()
            @assert $(map(esc,exs)...)
        end
    end
end

end
