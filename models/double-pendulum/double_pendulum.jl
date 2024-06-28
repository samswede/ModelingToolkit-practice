using ModelingToolkit

@connector PendulumPort begin
    # Across variables
    x(t)
    y(t)

    # Through variables (flow variables)
    T(t), [connect = Flow]

end

@mtkmodel Pendulum begin
    @parameters begin
        g=9.8
        L=1.0
    end
    @variables begin
        x(t)=0.0
        y(t)=0.0
        dx(t)=0.0
        dy(t)=0.0
        T(t)=0.0
    end
    @components begin
        port = PendulumPort()
    end
    @equations begin
        # connectors
        port.x ~ x
        port.y ~ y

        # physics
        D(x) ~ dx
        D(y) ~ dy
        D(dx) ~ T * x
        D(dy) ~ T * y - g
        x^2 + y^2 ~ L^2
    end
end
