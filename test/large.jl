Isotopologues("C494H776O148N136S4"; abtype = :total)
ti = time()
itl = Isotopologues("C494H776O148N136S4"; abtype = :total)
te = time()
if te - ti > 1
    @info "`Isotopologues` takes too much time for chemical ≈ 10 kDa. Skip tests for larger chemicals and transitions."
else
    ti = time()
    itl = Isotopologues("C1482H2328O444N408S12"; abtype = :total)
    te = time()
    if te - ti > 2
        @info "`Isotopologues` takes too much time for chemical ≈ 30 kDa. Skip tests for larger chemicals and transitions."
    else
        ti = time()
        itl = Isotopologues("C2964H4656O888N816S24"; abtype = :total)
        te = time()
        if te - ti > 5
            @info "`Isotopologues` takes too much time for chemical ≈ 50 kDa. Skip tests for larger chemicals and transitions."
        else
            itl = Isotopologues("C4940H7760O1480N1360S40"; abtype = :total)
            Isotopologues("C494H776O148N136S4" => "C247H388O74N68S2"; abtype = :total)
            ti = time()
            itll = Isotopologues("C494H776O148N136S4" => "C247H388O74N68S2"; abtype = :max)
            te = time()
            if te - ti > 2
                @info "`Isotopologues` takes too much time for transitions ≈ 10 kDa -> 5 kDa. Skip tests for larger transitions."
            else
                itll = TandemIsotopologues("C988H1552O296N272S8" => "C494H776O148N136S4"; abtype = :max)
                itll = Isotopologues("C988H1552O296N272S8" => "C494H776O148N136S4"; abtype = :max)
            end
        end
    end
end