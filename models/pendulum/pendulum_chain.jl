using ModelingToolkit

@mtkmodel ImmovablePin begin
    @parameters begin
        a=0.0
        b=0.0
    end
    @variables begin
        x(t)
        y(t)
        T(t)
    end
    @components begin
        port = PendulumPort()
    end
    @equations begin
        # Position of the immovable pin
        x ~ a
        y ~ b

        # Connector
        port.x ~ x
        port.y ~ y
        port.T ~ -T # The force is opposite to the tension in the pendulum

    end
end


@connector PendulumPort begin
    # Across variables
    x(t)
    y(t) 

    # Through variables (flow variables)
    T(t), [connect = Flow] # Tension in the pendulum

end

@connector MechanicalPort begin
    # Across variables
    dx(t) # Velocity in the x direction
    dy(t) # Velocity in the y direction

    # Through variables (flow variables)
    Fx(t), [connect = Flow] # Force in the x direction
    Fy(t), [connect = Flow] # Force in the y direction

end

@mtkmodel Pendulum begin
    @parameters begin
        g=9.8
        L=1.0
    end
    @variables begin
        x(t)
        y(t)
        dx(t)
        dy(t)
        T(t)
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
