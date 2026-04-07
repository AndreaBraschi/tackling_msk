function itemNames = getItemNames(componentSet)
    num_items = componentSet.getSize();
    itemNames = cell(1, num_items);

    for i = 1:num_items
        item = componentSet.get(i - 1);
        itemName = item.getName();
        itemNames{i} = char(itemName); % collect current itemName
    end
end
